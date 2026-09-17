function [xopt, fopt, exitflag, output] = optimize(funcname, x0, A, b, ...
    Aeq, beq, lb, ub, opt_struct, gradients)

% Written by Andrew Ning.  Feb 2016.
% FLOW Lab, Brigham Young University.

    % set options
    options = optimoptions('fmincon');
    names = fieldnames(opt_struct);
    for i = 1:length(names)
        options = optimoptions(options, names{i}, opt_struct.(names{i}));
    end

    % check if gradients provided
    if gradients
        options = optimoptions(options, 'GradObj', 'on', 'GradConstr', 'on');
    end

    % shared variables
    xcheck = 2*ub;
    cin = [];
    gcin = [];

    % run fmincon
    [xopt, fopt, exitflag, output] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);

    % ----- conversion functions from numpy to matlab ---------
    function [matrix] = matrixfromnumpy(array)
        data = double(py.array.array('d', py.numpy.nditer(array, pyargs('order', 'F'))));  % Add order='F' to get data in column-major order (as in Fortran 'F' and Matlab)
        data_size = cell2mat(cell(array.shape));
        matrix = reshape(data, data_size);  % No need for transpose, since we're retrieving the data in column major order
    end

    function [vector] = vectorfromnumpy(array)
        vector = double(py.array.array('d', array));
    end


    % ---------- Update objectives and constraints ------------------
    function [J, cin, gJ, gcin] = fupdate(x)
        eval(['output = py.', funcname, '(x);'])  % output is a cell array with {J, cin}
        xcheck = x;

        J = output{1};
        cin = vectorfromnumpy(output{2});

        if gradients
            gJ = vectorfromnumpy(output{3});
            gcin = matrixfromnumpy(output{4});
        else
            gJ = [];
            gcin = [];
        end
    end


    % ---------- Objective Function ------------------
    function [J, gJ] = obj(x)
        [J, cin, gJ, gcin] = fupdate(x);
    end
    % -------------------------------------------------

    % ------------- Constraints ------------------------
    function [c, ceq, gc, geq] = con(x)
        if any(x ~= xcheck)
            [J, cin, gJ, gcin] = fupdate(x);
        end
        c = cin;
        ceq = [];
        gc = gcin;
        geq = [];
    end
    % ------------------------------------------------


end
