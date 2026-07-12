classdef KeyedDedup < handle
%KEYEDDEDUP  O(1)-amortized solution dedup via a bucket hash. Match norm selectable.
%
%  Replaces the O(N) linear distance scan with a bucketed nearest-match lookup.
%  MATCH criterion (matchnorm):
%     'inf' :  max_k |sig_k(stored) - sig_k(cand)| < tol      (scale-invariant; matches the
%              per-variable Newton convergence tolerance -- recommended with vl = 1:numofvar)
%     'L1'  :  sum_k |sig_k(stored) - sig_k(cand)| < tol       (the classic sum; use only with
%              a SMALL sig, e.g. vl = a few indices, so the sum stays bounded as V grows)
%  sig = x(vl), optionally magnitude-folded (abs). In exact mode (useAbs=false) the antipode
%  (stored + cand) is also tested. With vl = 1:numofvar the match compares the full vectors
%  and never merges distinct solutions.
%
%  BUCKET KEY = sum over a small subset klist of the signature: key(s) = sum(s(klist)).
%  It is 1-Lipschitz under L1 and K-Lipschitz (K=numel(klist)) under the inf norm, so a
%  match lives within +/-1 bucket when the bucket width is tol ('L1') or K*tol ('inf').
%  A few variables keep their spread, so the key discriminates and buckets stay tiny.
%  The key affects ONLY speed; correctness is set entirely by the match.
%
%  Usage:
%    dd = KeyedDedup(1:numofvar, 1e-7, false, cap, klist, 'inf');  % full inf-norm match
%    dd = KeyedDedup(klist,      4e-7, false, cap, [],    'L1');   % few-var L1 match
%    dd.seed(Xcols);  [idx,isNew] = dd.query_or_add(z);

    properties
        vl; tol; useAbs; klist; bw; isL1;
        sig; bkey; n; cap;
        map;     % containers.Map  int64 bucket -> row vector of stored indices
    end
    methods
        function o = KeyedDedup(vl, tol, useAbs, cap0, klist, matchnorm)
            if nargin<4 || isempty(cap0); cap0 = 1024; end
            o.vl=vl(:).'; o.tol=tol; o.useAbs=logical(useAbs);
            if nargin<5 || isempty(klist); o.klist = 1:numel(o.vl); else; o.klist = klist(:).'; end
            if nargin<6 || isempty(matchnorm); matchnorm='inf'; end
            o.isL1 = strcmpi(matchnorm,'L1');
            if o.isL1; o.bw = o.tol; else; o.bw = numel(o.klist)*o.tol; end
            o.cap=max(cap0,16); o.sig=zeros(numel(o.vl),o.cap); o.bkey=zeros(1,o.cap);
            o.n=0; o.map=containers.Map('KeyType','int64','ValueType','any');
        end

        function seed(o, X)
            for j=1:size(X,2)
                s = X(:,j); s = s(o.vl); s = s(:); if o.useAbs; s = abs(s); end
                ks = sum(s(o.klist));
                o.n = o.n+1;
                if o.n>o.cap; o.cap=2*o.cap; o.sig(:,o.cap)=0; o.bkey(o.cap)=0; end
                o.sig(:,o.n)=s; o.bkey(o.n)=ks;
                b = int64(floor(ks/o.bw));
                if isKey(o.map,b); o.map(b)=[o.map(b) o.n]; else; o.map(b)=o.n; end
            end
        end

        function [idx,isNew] = query_or_add(o, z)
            s = z(o.vl); s = s(:);
            if o.useAbs; s = abs(s); end
            ks = sum(s(o.klist));
            idx = o.find_match(s, ks);
            if idx>0; isNew=false; return; end
            o.n = o.n+1;
            if o.n>o.cap; o.cap=2*o.cap; o.sig(:,o.cap)=0; o.bkey(o.cap)=0; end
            o.sig(:,o.n)=s; o.bkey(o.n)=ks;
            b = int64(floor(ks/o.bw));
            if isKey(o.map,b); o.map(b)=[o.map(b) o.n]; else; o.map(b)=o.n; end
            idx=o.n; isNew=true;
        end
    end
    methods (Access=private)
        function idx = find_match(o, s, ks)
            idx=0; w=o.tol;
            b0=int64(floor(ks/o.bw)); blist=[b0-1 b0 b0+1];
            if ~o.useAbs; bn=int64(floor(-ks/o.bw)); blist=[blist bn-1 bn bn+1]; end
            blist=unique(blist);
            for b=blist
                if isKey(o.map,b)
                    cols=o.map(b);
                    for j=cols
                        sj=o.sig(:,j);
                        if o.isL1
                            if sum(abs(sj-s))<w; idx=j; return; end
                            if ~o.useAbs && sum(abs(sj+s))<w; idx=j; return; end
                        else
                            if max(abs(sj-s))<w; idx=j; return; end
                            if ~o.useAbs && max(abs(sj+s))<w; idx=j; return; end
                        end
                    end
                end
            end
        end
    end
end
