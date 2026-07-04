function mpc = case7Salam
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
	1	1	-90     -30     0	0	1	1   	0	230	1	1.1	0.9;
	2	1	47.8     3.9	0	0	1	1       0	230	1	1.1	0.9;
	3	1	94.2    19.0	0	0	1	1       0	230	1	1.1	0.9;
    4	1	13.5    5.8     0	0	1	1       0	230	1	1.1	0.9;
    5	1	18.3    12.7    0	0	1	1   	0	230	1	1.1	0.9;
    6	1	 7.6     1.6    0	0	1	1   	0	230	1	1.1	0.9;
    7	3	0       0       0	0	1	1   	0	230	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	7	80	0	100	-100	1       100	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.082	0.192	0.0	250	250	250	0	0	1	-360	360;
	2	3	0.067	0.171	0.0	250	250	250	0	0	1	-360	360;
    2	5	0.058	0.176	0.0	250	250	250	0	0	1	-360	360;
    2	6	0.013	0.042	0.0	250	250	250	0	0	1	-360	360;
    3	4	0.024	0.100	0.0	250	250	250	0	0	1	-360	360;
    4	5	0.024	0.100	0.0	250	250	250	0	0	1	-360	360;
    5	6	0.057	0.174	0.0	250	250	250	0	0	1	-360	360;
    5	7	0.019	0.059	0.0	250	250	250	0	0	1	-360	360;
    6	7	0.054	0.223	0.0	250	250	250	0	0	1	-360	360;
];
