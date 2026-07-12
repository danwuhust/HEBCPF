function mpc = case5loop
%CASE4GS  Power flow data for 5 bus system loop system
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
	2	2	20  15	0	0	1	1	0	230	1	1.1	0.9;
	3	2	20  40	0	0	1	1	0	230	1	1.1	0.9;
    4	2	35   5	0	0	1	1	0	230	1	1.1	0.9;
    5	2	70  55	0	0	1	1	0	230	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	30	0   100	-100	1	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	25	0	100	-100	1	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
    3	55	0	100	-100	1	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
    4   50	0	100	-100	1	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
    5   45	0	100	-100	1	100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.01	0.10	0.0	250	250	250	0	0	1	-360	360;
	2	3	0.20	2.00	0.0	250	250	250	0	0	1	-360	360;
    3	4	0.20	2.00	0.0	250	250	250	0	0	1	-360	360;
    4	5	0.01	0.10	0.0	250	250	250	0	0	1	-360	360;
    1	5	0.01	0.10	0.0	250	250	250	0	0	1	-360	360;
%     3	4	0.02	0.10	0.0	250	250	250	0	0	1	-360	360;
];
