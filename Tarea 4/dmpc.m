function [accion_control, x_pred] = dmpc(dg_index, mediciones, parametros, pred_dg_vec, a_adj, parametros_mpc)
    % lambda1 penalización seguimiento de voltaje
    % lambda2 penalización control secundario
    % lambda3 penalización despacho económico
    % Ny horizonte predicción 
    % Nu horizonte acción de control
    % Mp coeficiente droop
    % Pmin_i restricción potencia mínima
    % Pmax_i restricción potencia máxima
    % Vmin_i restricción voltaje mínimo
    % Vmax_i restricción voltaje máximo


    dim_x = 4;
    dim_u = 1;

    %yalmip('clear')
    R = parametros(1);
    Mp = parametros(2);
    a_i = parametros(3);
    b_i = parametros(4);
    Pmax_i = parametros(5);
    V0_dc = parametros(6);
    Vmin_dc = parametros(7);
    Vmax_dc = parametros(8);

    Ny = parametros_mpc(1);
    Nu = parametros_mpc(2);
    lambda1 = parametros_mpc(3);
    lambda2 = parametros_mpc(4);
    lambda3 = parametros_mpc(5); 

    pred_size = 4*Ny;
    x_pred = zeros(1, pred_size);
    coder.varsize('x_pred', [1 pred_size])

    V0 = mediciones(1);
    i_medido = mediciones(2);
    P0 = mediciones(3);

    Vr0 = V0 - i_medido*R;

    Vmax_i = Vmax_dc;

    Vmin_i = Vmin_dc;
    Pmin_i = Vmin_i*Mp + V0_dc;
   
    aj = a_adj(dg_index, :);
    G = 1/R;
    J = 0;
    prom_con = sum(aj); 
    largo_vector_estado = Ny*dim_x+Nu*dim_u;
    H = zeros(largo_vector_estado);
    F = zeros(1, Ny*dim_x + Nu*dim_u);

    for i = 1:Ny
        H(i,i) = lambda1;
        F(i) = -2*V0;
        sum_despacho = 0;
        for j = 1:length(aj)
            if dg_index == j
                continue
            end
            preds_j = pred_dg_vec(:, j);
            sum_despacho = sum_despacho + aj(j)*preds_j(3*Ny + i);
        end
        H(Ny*3+i, Ny*3+i) = lambda3*(prom_con);
        F(Ny*3+1) = sum_despacho;

    end
    for i = 1:Nu
        H(4*Ny+i, 4*Ny+i) = lambda2;
    end

    % Restricciones de desigualdad
    lb = zeros(1, largo_vector_estado);
    lb(Ny : 2*Ny-1) = Vmin_i;
    lb(2*Ny : 3*Ny-1) = Pmin_i;

    ub = zeros(1, largo_vector_estado);
    ub(Ny : 2*Ny-1) = Vmax_i;
    ub(2*Ny : 3*Ny-1) = Pmax_i;

    C0 = G*(2*V0 - Vr0);
    n = 4*Ny + Nu;


    A = [];
    b = [];

    ny_1_vector = ones(1, Ny);
    % Dinámica voltaje
    v_mat = diag(ny_1_vector, Ny+1);
    v_mat = pad_mat(v_mat, n, Ny);

    v_mat_prev = diag(-ny_1_vector, Ny);
    v_mat_prev = pad_mat(v_mat_prev, n, Ny);

    p_mat = diag(-Mp*ny_1_vector, 2*Ny+1);
    p_mat = pad_mat(p_mat, n, Ny);

    p_mat_prev = Mp*diag(ny_1_vector, 2*Ny);
    p_mat_prev = pad_mat(p_mat_prev, n, Ny);

    u_prev = diag(-ny_1_vector, dim_x*Ny);
    u_prev = pad_mat(u_prev, n, Ny);

    restricciones_voltaje = v_mat + v_mat_prev + p_mat + p_mat_prev + u_prev;

    A = [A; restricciones_voltaje];
    b = [b; zeros(Ny, 1)];

    % Dinámica potencia
    p_mat = diag(ny_1_vector, 2*Ny+1);
    p_mat = pad_mat(p_mat, n, Ny);

    v_mat = diag(-C0*ny_1_vector, Ny+1);
    v_mat = pad_mat(v_mat, n, Ny);

    restricciones_potencia = p_mat + v_mat;

    A=[A; restricciones_potencia];
    b = [b; P0-C0*V0*ones(Ny, 1)];

    v_avg_mat = diag(prom_con*ny_1_vector);
    v_avg_mat = pad_mat(v_avg_mat, n, Ny);

    v_mat = diag(-ny_1_vector, Ny);
    v_mat = pad_mat(v_mat, n, Ny);

    restricciones_voltaje_promedio = v_avg_mat + v_mat;
    A = [A; restricciones_voltaje_promedio];
    %Restricciones estimación local de voltaje promedio

    for i = 1:Ny
        sum_v_pred = 0;
        for j = 1:length(aj)
            if dg_index == j
                continue
            end
            preds_j = pred_dg_vec(:, j);
            sum_v_pred = sum_v_pred + aj(j)*preds_j(Ny + i);
        end

        b = [b; sum_v_pred];
    end

    despacho_mat = diag(ny_1_vector, 3*Ny+1);
    despacho_mat = pad_mat(despacho_mat, n, Ny);

    p_mat = diag(-a_i*ny_1_vector, 2*Ny+2);
    p_mat = pad_mat(p_mat, n, Ny);

    restricciones_despacho = despacho_mat + p_mat;

    A = [A; restricciones_despacho];

    b = [b; b_i*ones(Ny, 1)];
    
    %%    
    x0 = zeros(length(H(:, 1)), 1);
    options = optimoptions('quadprog', 'Algorithm', 'active-set');
    [x, fval] = quadprog(H, F, [], [], A, b, lb, ub, x0, options);
    x = x';
    x_pred = x_pred + x(1:4*Ny);
    accion_control = x(dim_x*Ny + 1);
    
