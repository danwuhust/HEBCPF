function mpc = case3
%CASE4GS  Power flow data for 3 bus loop system - lossless, equal line
%parameters, no power injections.
%   Please see CASEFORMAT for details on the case file format.
%


%   MATPOWER
%   $Id: case3PV.m,v 1.0 2011/05/19 BCL $

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	 0	0	0	0	1	1	0	230	1	1.1	0.9;
	2	1	0  -300	0	0	1	1	120	230	1	1.1	0.9;
	3	1	0  -300	0	0	1	1	240	230	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	500	-500	1	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	300	100	-100	1	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     3   0	300	100	-100	1	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.0	1.00	0.0	2500	2500	2500	0	0	1	-360	360;
	2	3	0.0	1.00	0.0	2500	2500	2500	0	0	1	-360	360;
    3	1	0.0	1.00	0.0	2500	2500	2500	0	0	1	-360	360;
];
