function mpc = case5Salam_mod3
% CASE5 Salam
% Modified by Dan Wu to create a PV curve that fails standard holomorphic
% embedding techniques for approaching the singular point. 
% Increase the active power demand at bus 2 (equivalently decrease the
% distributed power generation at bus 2)


%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	20     10     0	0	1	1.0 	0	230	1	1.1	0.9;
	2	1	45     15	0	0	1	1       0	230	1	1.1	0.9;
	3	1	40      5	0	0	1	1       0	230	1	1.1	0.9;
    4	1	60     10     0	0	1	1       0	230	1	1.1	0.9;
    5	3	0       0       0	0	1	1.06	0	230	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	40	0	100	-100	1       100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
    5   0	0	100	-100	1.06	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.06	0.18	0.0	250	250	250	0	0	1	-360	360;
	1	3	0.06	0.18	0.0	250	250	250	0	0	1	-360	360;
    1	4	0.04	0.12	0.0	250	250	250	0	0	1	-360	360;
    1	5	0.02	0.06	0.0	250	250	250	0	0	1	-360	360;
    2	3	0.01	0.03	0.0	250	250	250	0	0	1	-360	360;
    2	5	0.08	0.24	0.0	250	250	250	0	0	1	-360	360;
    3	4	0.08	0.24	0.0	250	250	250	0	0	1	-360	360;
];
%% modification
% mpc.bus(2,4)=mpc.bus(2,4)+100; % decrease reactive power demand
mpc.bus(3,4)=mpc.bus(3,4)+80;
mpc.bus(4,4)=mpc.bus(4,4)-55;
% mpc.bus(2,3)=mpc.bus(2,3)+350;
mpc.bus(2,3)=mpc.bus(2,3)-550; % decrease distributed generation
mpc.bus(3,3)=mpc.bus(4,3)+200;
mpc.bus(4,3)=mpc.bus(4,3)+200;
