function build_mex()
%BUILD_MEX  Rebuild the two beneficial MEX functions for this solver.
%
% Produces:  Pade_Apprxmt_mex.mexw64   holomorphic_cont_tri_mex.mexw64
%
% REQUIREMENTS: MATLAB Coder + a configured C compiler (mex -setup C).
%
% WINDOWS BUILD GOTCHA (important):
%   If codegen fails with  "'<name>_mex.bat' is not recognized as an internal
%   or external command", the environment variable
%        NoDefaultCurrentDirectoryInExePath = 1
%   is set, which stops cmd.exe from running codegen's generated .bat from its
%   build folder. Fix: unset it BEFORE launching MATLAB, e.g. from a shell
%        set "NoDefaultCurrentDirectoryInExePath="   &&  matlab ...
%   (bash:  unset NoDefaultCurrentDirectoryInExePath; matlab ... )
%   The compiler itself is fine; only codegen's build step is affected.
%
% The types below treat 'degree' as a compile-time constant (act=15). If you
% change degree.act / numerator / denominator in main.m, regenerate examples.mat
% (any quick run that saves the degree/I/Start structs) and rebuild.

addpath(pwd);
S  = load('examples.mat');
DC = coder.Constant(S.degree);

% I.pv / I.pq may be Nx1 or empty 0x0 (e.g. case3 has no PV bus) -> 2-D variable
Itype = coder.typeof(S.I);
Itype.Fields.pv = coder.typeof(0,[Inf Inf],[1 1]);
Itype.Fields.pq = coder.typeof(0,[Inf Inf],[1 1]);

Stype = coder.typeof(S.Start);
Stype.Fields.V0 = coder.typeof(1i,[Inf 1],[1 0]);
Stype.Fields.W0 = coder.typeof(1i,[Inf 1],[1 0]);
Stype.Fields.P0 = coder.typeof(0 ,[Inf 1],[1 0]);
Stype.Fields.Q0 = coder.typeof(0 ,[Inf 1],[1 0]);

vecR = coder.typeof(0 ,[Inf 1]  ,[1 0]);
matR = coder.typeof(0 ,[Inf Inf],[1 1]);
sc   = 0;

codegen Pade_Apprxmt         -args {matR,matR,DC,sc}                                  -o Pade_Apprxmt_mex
codegen holomorphic_cont_tri -args {vecR,sc,vecR,vecR,matR,matR,Stype,Itype,DC,sc}    -o holomorphic_cont_tri_mex
fprintf('Built Pade_Apprxmt_mex and holomorphic_cont_tri_mex.\n');
end
