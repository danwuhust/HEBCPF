function [Ma,ra,I] = quadr_matrix(Mp,Mq,mpc,bus_n,gen_n,numofvar)
%QUADR_MATRIX Assemble sparse PV, PQ and reference-bus equations.

[PQ,PV,REF,~,BUS_I,BUS_TYPE,PD,QD] = idx_bus;
[GEN_BUS,PG,~,~,~,VG,~,GEN_STATUS] = idx_gen;

bus = mpc.bus;
gen = mpc.gen;
if size(bus,1) ~= bus_n || size(gen,1) ~= gen_n
    error('HEBCPF:quadr_matrix:DimensionMismatch', ...
        'bus_n and gen_n must match the MATPOWER case.');
end
if numofvar ~= 2*bus_n-1
    error('HEBCPF:quadr_matrix:VariableCount', ...
        'numofvar must equal 2*bus_n-1.');
end
if ~isequal(bus(:,BUS_I),(1:bus_n)')
    error('HEBCPF:quadr_matrix:ExternalNumbering', ...
        'Use ext2int before constructing HEBCPF matrices.');
end
if numel(Mp) ~= bus_n || numel(Mq) ~= bus_n
    error('HEBCPF:quadr_matrix:MatrixCount', ...
        'Mp and Mq must contain one matrix per bus.');
end

types = bus(:,BUS_TYPE);
unsupported = find(~ismember(types,[PQ PV REF]),1);
if ~isempty(unsupported)
    error('HEBCPF:quadr_matrix:UnsupportedBusType', ...
        'Bus %d has unsupported MATPOWER bus type %d.', ...
        bus(unsupported,BUS_I),types(unsupported));
end

I.pv = find(types == PV);
I.pq = find(types == PQ);
I.slack = find(types == REF);
if numel(I.slack) ~= 1
    error('HEBCPF:quadr_matrix:ReferenceBusCount', ...
        'HEBCPF requires exactly one reference bus; found %d.',numel(I.slack));
end
I.slack = I.slack(1);
I.pv = I.pv(:);
I.pq = I.pq(:);

npg = numel(I.pv);
npd = numel(I.pq);
nq = npd;
nv = npg;
Ma = cell(numofvar,1);
ra = zeros(numofvar,1);
ref_imag = bus_n + I.slack;
eq = 0;

eq = eq + 1;
Ma{eq} = voltage_matrix(I.slack,bus_n,ref_imag);
ra(eq) = voltage_target(gen,bus(I.slack,BUS_I),GEN_BUS,VG,GEN_STATUS);

for k = 1:nv
    i = I.pv(k);
    eq = eq + 1;
    Ma{eq} = voltage_matrix(i,bus_n,ref_imag);
    ra(eq) = voltage_target(gen,bus(i,BUS_I),GEN_BUS,VG,GEN_STATUS);
end

for k = 1:nq
    i = I.pq(k);
    eq = eq + 1;
    Ma{eq} = remove_reference_imag(Mq{i},ref_imag);
    ra(eq) = -bus(i,QD)/mpc.baseMVA;
end

for k = 1:npg
    i = I.pv(k);
    rows = online_generators(gen,bus(i,BUS_I),GEN_BUS,GEN_STATUS);
    eq = eq + 1;
    Ma{eq} = remove_reference_imag(Mp{i},ref_imag);
    ra(eq) = (sum(gen(rows,PG))-bus(i,PD))/mpc.baseMVA;
end

for k = 1:npd
    i = I.pq(k);
    eq = eq + 1;
    Ma{eq} = remove_reference_imag(Mp{i},ref_imag);
    ra(eq) = -bus(i,PD)/mpc.baseMVA;
end

if eq ~= numofvar
    error('HEBCPF:quadr_matrix:EquationCount', ...
        'Constructed %d equations, expected %d.',eq,numofvar);
end
end

function A = voltage_matrix(bus_idx,bus_n,ref_imag)
A = sparse([bus_idx; bus_n+bus_idx],[bus_idx; bus_n+bus_idx],1,2*bus_n,2*bus_n);
A = remove_reference_imag(A,ref_imag);
end

function A = remove_reference_imag(A,ref_imag)
keep = true(size(A,1),1);
keep(ref_imag) = false;
A = sparse(A(keep,keep));
end

function rows = online_generators(gen,bus_id,GEN_BUS,GEN_STATUS)
rows = find(gen(:,GEN_BUS) == bus_id & gen(:,GEN_STATUS) > 0);
if isempty(rows)
    error('HEBCPF:quadr_matrix:MissingGenerator', ...
        'PV/reference bus %d has no online generator.',bus_id);
end
end

function target = voltage_target(gen,bus_id,GEN_BUS,VG,GEN_STATUS)
rows = online_generators(gen,bus_id,GEN_BUS,GEN_STATUS);
setpoints = gen(rows,VG);
tolerance = 1e-10*max(1,abs(setpoints(1)));
if any(abs(setpoints-setpoints(1)) > tolerance)
    error('HEBCPF:quadr_matrix:InconsistentVoltageSetpoints', ...
        'Online generators at bus %d have inconsistent VG setpoints.',bus_id);
end
target = setpoints(1)^2;
end
