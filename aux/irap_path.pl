% =========================================================
% Copyright 2012-2013,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
%
% This file is part of iRAP.
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with iRAP.  If not, see <http://www.gnu.org/licenses/>.
%
%
%    $Id: scripts/irap Nuno Fonseca Sat Jan 5 20:56:39 2013$
% =========================================================


:-prolog_initialization(go).

% checks if a path (Mapper, Quant, Norm, DE) is valid
% or
% produces a plot with all the possible paths

go:-
	unix(argv(Args)),
	handle_args(Args),
	halt.

handle_args([File]):-!,
    start_graph(File).

handle_args([M,Q,N,D]):-!,
    (valid_combination([M,Q,N,D])->
     format("valid~n",[])
    ;
     format("invalid~n",[])
    ).
handle_args(_):-
    format("ERROR! usage: irap_paths [ FILENAME | Mapper Quant Norm DE ]~n",[]).

valid_combination([Map,QR,QN,DE]):-
    m(Map,_,_),
    qr(QR,m(Map),_),
    qn(QN,qr(QR),_),
    de(DE,(qr(QR),qn(QN)),_).

member(X,[X|_]).
member(X,[_|R]):-
    member(X,R).

% generate a graphviz file with all possible combinations
start_graph(File):-
	format('Generating graphviz file ~w~n',[File]),
	open(File,write,S),!,
	eraseall(r_graph_fd),
	recorda(r_graph_fd,S,_),
        graph_format('digraph { outputorder=edgesfirst;~n',[]),
        graph_format('overlap=scale;~n~n',[]),
%        graph_format('splines=ortho;~n~n',[]),
        graph_format('splines=polyline;~n~n',[]),
        graph_format('model=circuit;~n',[]),
%        graph_format('model=subset;~n',[]),
	save_nodes,
%	graph_format('cluster_m -> cluster_qr [style="invis"]~n',[]),
%	graph_format('cluster_qr -> cluster_qn [style="invis"]~n',[]),
%	graph_format('cluster_qn -> cluster_de [style="invis"]~n',[]),		
	save_edges,
	end_graph.
get_fd(FD):-
	recorded(r_graph_fd,FD,_).

end_graph:-
        graph_format('}~n',[]).

graph_format(A,B):-
	get_fd(FD),
        format(FD,A,B),
        !.

start_cluster(Type):-
        type2label(Type,Label,Style),
	atom_concat(['cluster_',Type],NTag),
	graph_format('subgraph ~w { ~n',[NTag]),
        graph_format('center=true;~n',[]),
	graph_format(' label=\"~w\";~n',[Label]),
	graph_format(' ~w~n',[Style]).

end_cluster:-
	graph_format('} ~n',[]).

type2label(m,'Mappers','style=filled; color=lightgrey;').
type2label(qr, 'Quantification','style=filled; color=lightgrey;').
type2label(qn,'Normalization','style=filled; color=lightgrey;').
type2label(de, 'DE','style=filled; color=lightgrey;').

%node(+nodeType,NodeName,Label,Style)
node(m,N,N,'[color=salmon2,style=filled]'):-
    m(N,_,_).
node(qr,N,N,'[color=coral3,style=filled]'):-
    qr(N,_,_).
node(qn,NDE,N,Style):-
    qn(N,_,_),
    atom_concat(['[color=cyan,style=filled,label=',N,']'],Style),
    atom_concat([N,qn],NDE).
node(de,NDE,N,Style):-
    de(N,_,_),
    atom_concat(['[color=lightsteelblue1,style=filled,label=',N,']'],Style),
    atom_concat([N,de],NDE).

save_nodes:-
    member(T,[m,qr,qn,de]),
    (start_cluster(T);end_cluster,fail),
    all(Label,node(T,_NodeName,Label,_),Ls),
    member(L,Ls),
    once(node(T,N,L,Attrs)),
    graph_format('~w ~w;~n',[N,Attrs]),
    fail.
save_nodes.

save_edges:-
    all((QR1,Map1),(m(Map1,_,_),qr(QR1,m(Map1),_)),L),
    member((QR,Map),L),
    graph_format('~w -> ~w;~n',[Map,QR]),
    fail.
save_edges:-
    all((QN1,QR1),(qr(QR1,_,_),qn(QN1,qr(QR1),_)),L),
    member((QN,QR),L),
    once(node(qn,QNE,QN,_)),
    graph_format('~w -> ~w;~n',[QR,QNE]),
    fail.

save_edges:-
    all((DE1,QN1),(qn(QN1,_,_),de(DE1,(_,qn(QN1)),_)),L),
    member((DE,QN),L),
    once(node(qn,QNE,QN,_)),
    once(node(de,NDE,DE,_)),
    graph_format('~w -> ~w;~n',[QNE,NDE]),
    fail.

save_edges:-
    fail,!,
    all((DE1,QR1),(qr(QR1,_,_),de(DE1,(qr(QR1),_),_)),L),
    member((DE,QR),L),
    once(node(qr,QRE,QR,_)),
    once(node(de,NDE,DE,_)),
    graph_format('~w -> ~w;~n',[QRE,NDE]),
    fail.
save_edges.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reference(genome).
%reference(transcriptome).
%junctions

m('tophat1',_,'').
m('tophat2',_,'').
m('smalt',_,'').
m('gsnap',_,'').
m('soapsplice',_,'').
m('bwa1',_,'').
m('bwa2',_,'').
m('bowtie1',_,'').
m('bowtie2',_,'').
m('gem',_,'').
m('gems',_,'').%RNA-version
m('star',_,'').
m('osa',_,'').

qr('htseq1',m(_),'Only requires the NH flag defined').
qr('htseq2',m(_),'Only requires the NH flag defined').
qr('basic',m(_),'').
qr('cufflinks1',m(M),'BAM flags...'):-m(M,_,_),not member(M,[soapsplice]).
qr('cufflinks1_nd',m(M),'BAM flags...'):-m(M,_,_),not member(M,[soapsplice]).
qr('cufflinks2',m(M),'BAM flags...'):-m(M,_,_),not member(M,[soapsplice]).
qr('cufflinks2_nd',m(M),'BAM flags...'):-m(M,_,_),not member(M,[soapsplice]).
qr('flux_cap',m(_),'').
qr('scripture',m(_),'').
qr('ireckon',m(_),'').
qr('bitseq',m(_),'').
qr('rsem',m(_),'').
qr('nurd',m(_),'').
qr('isoem',m(_),'').
qr('sailfish',m(_),'').

qr(none,m(_),'').
%qr('mmseq',m(_),'').


qn(cufflinks1,qr(cufflinks1),_).
qn(cufflinks2,qr(cufflinks2),_).
qn(cufflinks1_nd,qr(cufflinks1_nd),_).
qn(cufflinks2_nd,qr(cufflinks2_nd),_).
qn(nurd,qr(nurd),_).
qn(flux_cap,qr(flux_cap),_).
qn(deseq,qr(QR),_):-member(QR,[flux_cap,basic,htseq1,htseq2]).
qn(none,qr(_),_).


de(deseq,(qr(QR),qn(QN)),_):-member(QR,[htseq1,htseq2,basic,flux_cap]),member(QN,[deseq,flux_cap,none]).
de(dexseq,(qr(QR),qn(QN)),_):-member(QR,[htseq1,htseq2,basic,flux_cap]),member(QN,[deseq,flux_cap,none]).
de(bayseq,(qr(QR),qn(QN)),_):-member(QR,[htseq1,htseq2,basic,flux_cap]),member(QN,[deseq,flux_cap,none]).
de(edger,(qr(QR),qn(QN)),_):-member(QR,[htseq1,htseq2,basic,flux_cap]),member(QN,[deseq,flux_cap,none]).
de(voom,(qr(QR),qn(QN)),_):-member(QR,[htseq1,htseq2,basic,flux_cap]),member(QN,[deseq,flux_cap,none]).
de(cuffdiff1,(qr(QR),qn(QR)),_):-member(QR,[cufflinks1,cufflinks2]).
de(cuffdiff2,(qr(QR),qn(QR)),_):-member(QR,[cufflinks1,cufflinks2]).
de(cuffdiff1_nd,(qr(QR),qn(QR)),_):-member(QR,[cufflinks1_nd,cufflinks2_nd]).
de(cuffdiff2_nd,(qr(QR),qn(QR)),_):-member(QR,[cufflinks1_nd,cufflinks2_nd]).
de(none,(qr(_),qn(_)),_).
%de(edger,(qr(QR),_),_):-member(QR,[htseq1,htseq2,basic,flux_cap]).



