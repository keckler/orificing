run('orificeInput.m');

fprintf('beginning\n');

%%%%%
%pre-optimization manipulation
%%%%%

fprintf('reading in power files\n');

Q = readQ(powerDetectorFiles);
lengthQ_original = length(Q);
[Q, map] = formMap(Q, assemblyPowerThreshold);
Q_ave = sum(Q,2)/length(powerDetectorFiles); %divide by number of steps to get average assembly power over cycle
alpha = Q/cp; %kg-K/s

fprintf('manipulating data\n');

nass = length(Q_ave); %number of assemblies in problem
nsteps = length(powerDetectorFiles);
npossflows = length(x); %number of possible flowrates specified as data
nvars = nass+nass*npossflows+npossflows; %number of total variables

%build tensor of discretized dT
Omega = zeros(nass, npossflows, nsteps);
for k = 1:nsteps
    Omega(:,:,k) = alpha(:,k)*(1./x);
end

%find which assemblies have adjacents and how many adjacent pairs
adjacentAssemblies = findAdjacentAssemblies(nass, map, lengthQ_original);
nadj = nnz(adjacentAssemblies); %number of adjacent assembly pairs

%find how many constraints
nineqs = nass*nsteps + 2*nsteps + nass + 2*nadj*nsteps + npossflows*nass; %number of inequality constraints
neqs = nass + nass; %number of equality constraints

fprintf('initializing problem structures\n');

%initialize constraint matrices and vectors
nelements_Aineq = nass+nsteps*nadj*npossflows*2+nass*npossflows+nsteps*nass+nsteps*nass*2+nass*npossflows;
Aineq_i = zeros(nelements_Aineq,1);
Aineq_j = zeros(nelements_Aineq,1);
Aineq_v = zeros(nelements_Aineq,1);
bineq = zeros(nineqs, 1);
nelements_Aeq = (npossflows+1)*nass + npossflows*nass;
Aeq_i = zeros(nelements_Aeq,1);
Aeq_j = zeros(nelements_Aeq,1);
Aeq_v = zeros(nelements_Aeq,1);
beq = zeros(neqs, 1);
c = zeros(nvars,1);

%specify which variables are continuous versus integer/binary
ctype = '';
for i = 1:nvars
    if i <= nass %mdot variables
        ctype(end+1) = 'C';
    elseif i <= nass+nass*npossflows %delta variables
        ctype(end+1) = 'B';
    else %beta variables
        ctype(end+1) = 'C';
    end
end

%%%%%
%specify the constraint matrices
%%%%%

fprintf('determining the constraints\n');

idx_Aineq = 1; %keep track of the number of elements in Aineq matrix
idx_Aeq = 1;

%max outlet temp (constraint 1)
for k = 1:nsteps
    for i = 1:nass
        Aineq_i(idx_Aineq) = (k-1)*nass+i;
        Aineq_j(idx_Aineq) = i;
        Aineq_v(idx_Aineq) = -1;
        idx_Aineq = idx_Aineq+1;
        
        bineq((k-1)*nass+i) = -alpha(i,k)/dT_max;
    end
end

%mixed outlet temp (constraint 2)
for k = 1:nsteps
    for i = 1:nass
        Aineq_i(idx_Aineq) = nsteps*nass+k;
        Aineq_j(idx_Aineq) = i;
        Aineq_v(idx_Aineq) = -1;
        idx_Aineq = idx_Aineq+1;
        
        Aineq_i(idx_Aineq) = nsteps*nass+nsteps+k;
        Aineq_j(idx_Aineq) = i;
        Aineq_v(idx_Aineq) = 1;
        idx_Aineq = idx_Aineq+1;
    end
    bineq(nsteps*nass+k) = -sum(alpha(:,k))/(T_out_bar+T_out_bar_tol-T_in);
    bineq(nsteps*nass+nsteps+k) = sum(alpha(:,k))/(T_out_bar-T_out_bar_tol-T_in);
end

%max flow (constraint 3)
for i = 1:nass
    Aineq_i(idx_Aineq) = nsteps*nass+2*nsteps+i;
    Aineq_j(idx_Aineq) = i;
    Aineq_v(idx_Aineq) = 1;
    idx_Aineq = idx_Aineq+1;
    
    bineq(nsteps*nass+2*nsteps+i) = v_max*rho*A_flow;
end

%adjacent outlet temp (constraint 4)
constraint_idx = 1;
for k = 1:nsteps
    for i = 1:nass
        for ip = adjacentAssemblies(i,:)
            if ip > 0
                Aineq_i(idx_Aineq:idx_Aineq+npossflows-1) = nsteps*nass+2*nsteps+nass+constraint_idx;
                Aineq_j(idx_Aineq:idx_Aineq+npossflows-1) = [nass+(i-1)*npossflows+1:nass+(i-1)*npossflows+npossflows];
                Aineq_v(idx_Aineq:idx_Aineq+npossflows-1) = Omega(i,:,k);
                Aineq_i(idx_Aineq+npossflows:idx_Aineq+2*npossflows-1) = nsteps*nass+2*nsteps+nass+constraint_idx;
                Aineq_j(idx_Aineq+npossflows:idx_Aineq+2*npossflows-1) = [nass+(ip-1)*npossflows+1:nass+(ip-1)*npossflows+npossflows];
                Aineq_v(idx_Aineq+npossflows:idx_Aineq+2*npossflows-1) = -Omega(ip,:,k);
                idx_Aineq = idx_Aineq+2*npossflows;
                
                Aineq_i(idx_Aineq:idx_Aineq+npossflows-1) = nsteps*nass+2*nsteps+nass+constraint_idx+1;
                Aineq_j(idx_Aineq:idx_Aineq+npossflows-1) = [nass+(i-1)*npossflows+1:nass+(i-1)*npossflows+npossflows];
                Aineq_v(idx_Aineq:idx_Aineq+npossflows-1) = -Omega(i,:,k);
                Aineq_i(idx_Aineq+npossflows:idx_Aineq+2*npossflows-1) = nsteps*nass+2*nsteps+nass+constraint_idx+1;
                Aineq_j(idx_Aineq+npossflows:idx_Aineq+2*npossflows-1) = [nass+(ip-1)*npossflows+1:nass+(ip-1)*npossflows+npossflows];
                Aineq_v(idx_Aineq+npossflows:idx_Aineq+2*npossflows-1) = Omega(ip,:,k);
                idx_Aineq = idx_Aineq+2*npossflows;
                
                constraint_idx = constraint_idx + 2;
            end
        end
    end
end
bineq(nsteps*nass+2*nsteps+nass+1:nsteps*nass+2*nsteps+nass+2*nadj*nsteps) = xi;

%number of groups (constraint 8)
for i = 1:nass
    for j = 1:npossflows
        Aineq_i(idx_Aineq) = nass*nsteps+2*nsteps+nass+2*nadj*nsteps+(i-1)*npossflows+j;
        Aineq_j(idx_Aineq) = nass+(i-1)*npossflows+j;
        Aineq_v(idx_Aineq) = 1;
        idx_Aineq = idx_Aineq+1;
        
        Aineq_i(idx_Aineq) = nass*nsteps+2*nsteps+nass+2*nadj*nsteps+(i-1)*npossflows+j;
        Aineq_j(idx_Aineq) = nass+nass*npossflows+j;
        Aineq_v(idx_Aineq) = -1;
        idx_Aineq = idx_Aineq+1;
    end
end
bineq(nsteps*nass+2*nsteps+nass+2*nadj*nsteps+1:end) = 0;

%select flow from discretized flows (constraint 5)
for i = 1:nass
    Aeq_i(idx_Aeq) = i;
    Aeq_j(idx_Aeq) = i;
    Aeq_v(idx_Aeq) = -1;
    idx_Aeq = idx_Aeq+1;
    for j = 1:npossflows
        Aeq_i(idx_Aeq) = i;
        Aeq_j(idx_Aeq) = nass+(i-1)*npossflows+j;
        Aeq_v(idx_Aeq) = x(j);
        idx_Aeq = idx_Aeq+1;
    end
    beq(i) = 0;
end

%select only one flowrate for a channel (constraint 6)
for i = 1:nass
    for j = 1:npossflows
        Aeq_i(idx_Aeq) = nass+i;
        Aeq_j(idx_Aeq) = nass+(i-1)*npossflows+j;
        Aeq_v(idx_Aeq) = 1;
        idx_Aeq = idx_Aeq+1;
    end
    beq(nass+i) = 1;
end

fprintf('building the constraint matrices\n');

%build sparse matrices for Aineq and Aeq
Aineq = sparse(Aineq_i, Aineq_j, Aineq_v, nineqs, nvars);
Aeq = sparse(Aeq_i, Aeq_j, Aeq_v, neqs, nvars);

%objective
for j = 1:npossflows
    c(nass+nass*npossflows+j) = 1;
end

%variable bounds
%lb = -Inf*ones(numAss*numBatches+6,1);
%ub = Inf*ones(numAss*numBatches+6,1);

%print general problem parameters
fprintf('number of constraints = %i\n', neqs+nineqs);
fprintf('             equality = %i\n', neqs);
fprintf('           inequality = %i\n', nineqs);
fprintf('number of variables = %i\n', nvars);
fprintf('            integer = %i\n', sum(ctype == 'I'));
fprintf('             binary = %i\n', sum(ctype == 'B'));
fprintf('         continuous = %i\n', sum(ctype == 'C'));

%%%%%
%solve
%%%%%

fprintf('passing to cplex for solve\n');

%opts=cplexoptimset('display','on');
opts=cplexoptimset('display','iter');
opts.exportmodel = 'cplex_output.lp';
[solutionvector, objval, status, output] = cplexmilp(c, Aineq, bineq, Aeq, beq, [], [], [], [], [], ctype,[],opts);
fprintf('exit status = %i\n', status);
fprintf('solution time = %f\n', output.time);

%%%%%
%post-process
%%%%%

m = solutionvector(1:nass);
Tout = T_in + alpha./m;

%check max outlet temp and mixed outlet temp constraints
for k = 1:nsteps
    if sum(Tout(:,k) > T_in+dT_max) > 0
        fprintf('error: an assembly in step %i violates the outlet temperature constraint\n', k);
    end
    
    if sum((T_in + sum(alpha(:,k))/sum(m)) > T_out_bar+T_out_bar_tol) > 0
        fprintf('error: mixed outlet temperature violates constraint in step %i\n', k);
    elseif sum((T_in + sum(alpha(:,k))/sum(m)) < T_out_bar-T_out_bar_tol) > 0
        fprintf('error: mixed outlet temperature violates constraint in step %i\n', k);
    end
end

%check maximum flow velocity constraint
if sum(m/rho/A_flow > v_max) > 0
    fprintf('error: an assembly violates the maximum velocity constraint\n');
end

%check adjacent outlet temp constraint
