(* ::Package:: *)

(* ::Input:: *)
(*(* I. calculate the field equations using spherical decomposition; 35 mins *)*)


(* ::Input:: *)
(*(* Cell 1: represent Subscript[Y, l,m+2] and Subscript[Y, l,m+3] using Subscript[Y, l,m], Subscript[Y, l,m+1], and Subscript[Y, l,m+2] *)*)


(* ::Input::Initialization:: *)
lsqylm=Simplify[D[D[SphericalHarmonicY[l,m,th,ph],th],th]+1/(Sin[th])^2*D[D[SphericalHarmonicY[l,m,th,ph],ph],ph]+Cot[th]*D[SphericalHarmonicY[l,m,th,ph],th]]/.{SphericalHarmonicY[l,m,th,ph]->Sqrt[Gamma[l-m+1]/Gamma[l+m+1]]*Exp[I*m*ph]*plm,SphericalHarmonicY[l,m+1,th,ph]->Sqrt[Gamma[l-m]/Gamma[l+m+2]]*Exp[I*(m+1)*ph]*plm1,SphericalHarmonicY[l,m+2,th,ph]->Sqrt[Gamma[l-m-1]/Gamma[l+m+3]]*Exp[I*(m+2)*ph]*plm2};lsqylm=Simplify[D[D[SphericalHarmonicY[l,m,th,ph],th],th]+1/(Sin[th])^2*D[D[SphericalHarmonicY[l,m,th,ph],ph],ph]+Cot[th]*D[SphericalHarmonicY[l,m,th,ph],th]]/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2};ylm2iter=Simplify[ylm2/.Solve[lsqylm==-l*(l+1)*ylm,ylm2][[1]],Assumptions->{l>m>0}];lpylm=Exp[I*ph]*(D[SphericalHarmonicY[l,m,th,ph],th]+I*Cot[th]*D[SphericalHarmonicY[l,m,th,ph],ph]);lsqlpylm=D[D[lpylm,th],th]+(1/(Sin[th])^2)*D[D[lpylm,ph],ph]+Cot[th]*D[lpylm,th]/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2,SphericalHarmonicY[l,m+3,th,ph]->ylm3};ylm3iter=ylm3/.Solve[lsqlpylm==-l*(l+1)*Sqrt[(l-m)(l+m+1)]*ylm1,ylm3][[1]];ylm3iter=Simplify[ylm3iter/.ylm2->ylm2iter,Assumptions->{l>m>0}];


(* ::Input:: *)
(*(* Cell 2: define the zeroth-order quantities *)*)


(* ::Input::Initialization:: *)
n=4;coord={t,r,th,ph};metric0={{-bf[r],0,0,0},{0,af[r],0,0},{0,0,r^2,0},{0,0,0,r^2*(Sin[th])^2}};inversemetric0=Inverse[metric0];chris0=Simplify[Table[(1/2)*Sum[inversemetric0[[i,s]]*(D[metric0[[s,k]],coord[[j]]]+D[metric0[[j,s]],coord[[k]]]-D[metric0[[j,k]],coord[[s]]]),{s,1,n}],{i,1,n},{j,1,n},{k,1,n}]];u4velup0=(1/Sqrt[-(metric0[[1,1]])])*{1,0,0,0};u4veldown0=Table[Sum[metric0[[i,j]]*u4velup0[[j]],{j,1,n}],{i,1,n}];


(* ::Input:: *)
(*u4velup0*)


(* ::Input:: *)
(*u4veldown0*)


(* ::Input:: *)
(*(* Cell 3: define the even perturbation functions: {eh0f, eh1f, eh2f, ekf, w0f, w1f, vf, debphf} \[LongLeftRightArrow] {Subscript[H, 0], Subscript[H, 1], Subscript[H, 2], K, Subscript[W, 0], Subscript[W, 1], V, X}; po denotes the perturbation parameter *)*)


(* ::Input::Initialization:: *)
eh0=eh0f[t,r];eh1=eh1f[t,r];eh2=eh2f[t,r];ek=ekf[t,r];emetric1=po*{{bf[r]*eh0*SphericalHarmonicY[l,m,th,ph],eh1*SphericalHarmonicY[l,m,th,ph],0,0},{eh1*SphericalHarmonicY[l,m,th,ph],af[r]*eh2*SphericalHarmonicY[l,m,th,ph],0,0},{0,0,r^2*ek*SphericalHarmonicY[l,m,th,ph],0},{0,0,0,r^2*(Sin[th])^2*ek*SphericalHarmonicY[l,m,th,ph]}};fluw0=w0f[t,r];fluw1=w1f[t,r];fluv=vf[t,r];edisxdown=po*{fluw0*SphericalHarmonicY[l,m,th,ph],fluw1*SphericalHarmonicY[l,m,th,ph],fluv*D[SphericalHarmonicY[l,m,th,ph],th],fluv*D[SphericalHarmonicY[l,m,th,ph],ph]};edisxup=Table[Sum[edisxdown[[i1]]*inversemetric0[[i,i1]],{i1,1,n}],{i,1,n}];edebph=po*debphf[t,r]*SphericalHarmonicY[l,m,th,ph];


(* ::Input:: *)
(*(* Cell 4: setup the metric and fluid quantities  *)*)


(* ::Input::Initialization:: *)
metricp=emetric1;metric=metric0+metricp;inversemetric=Simplify[Normal[Series[Inverse[metric],{po,0,1}]],Assumptions->l>m>0];disxup=edisxup;disxdown=edisxdown;debph=edebph;bph=bphf[r]+debph;debph=edebph;bph=bphf[r]+debph;


(* ::Input:: *)
(*(* Lie derivatives *)*)


(* ::Input::Initialization:: *)
liedisxup=Simplify[Table[Sum[u4velup0[[s1]]*D[disxup[[i]],coord[[s1]]]-disxup[[s1]]*D[u4velup0[[i]],coord[[s1]]],{s1,1,n}]+Sum[chris0[[i,s1,s2]]*(u4velup0[[s1]]*disxup[[s2]]-u4velup0[[s2]]*disxup[[s1]]),{s1,1,n},{s2,1,n}],{i,1,n}],Assumptions->l>m>0];u4velupp=Simplify[Table[liedisxup[[i]]+u4velup0[[i]]*(Sum[liedisxup[[s1]]*u4veldown0[[s1]],{s1,1,n}]+1/2*Sum[u4velup0[[s1]]*u4velup0[[s2]]*metricp[[s1,s2]],{s1,1,n},{s2,1,n}]),{i,1,n}],Assumptions->l>m>0];u4velup=u4velup0+u4velupp;u4veldown=Simplify[Normal[Series[Table[Sum[u4velup[[s1]]*metric[[s1,i]],{s1,1,n}],{i,1,n}],{po,0,1}]],Assumptions->l>m>0];


(* ::Input:: *)
(*u4velup*)


(* ::Input:: *)
(*(* fluid perturbations \[DifferenceDelta]V, \[Delta]\[Epsilon], \[Delta]p *)*)


(* ::Input::Initialization:: *)
devv=Simplify[1/2*(1/af[r]*metricp[[2,2]]+1/r^2*metricp[[3,3]]+1/(r^2*(Sin[th])^2)*metricp[[4,4]])+Sum[D[disxup[[s1]],coord[[s1]]],{s1,2,n}]+1/2*Sum[inversemetric0[[s1,s2]]*D[metric0[[s1,s2]],coord[[s3]]]*disxup[[s3]],{s1,2,n},{s2,2,n},{s3,2,n}],Assumptions->l>m>0];deep=-(epf[r]+pf[r])*devv-disxup[[2]]*epf'[r];dep=-gaf[r]*pf[r]*devv-disxup[[2]]*pf'[r];


(* ::Input:: *)
(*dep*)


(* ::Input:: *)
(*(* calculate quantites in the field equation *)*)


(* ::Input:: *)
(*(* Ricci tensor *)*)


(* ::Input::Initialization:: *)
chris=Table[(1/2)*Sum[inversemetric[[i,s]]*(D[metric[[s,k]],coord[[j]]]+D[metric[[j,s]],coord[[k]]]-D[metric[[j,k]],coord[[s]]]),{s,1,n}],{i,1,n},{j,1,n},{k,1,n}];riemann=Table[D[chris[[i,j,ll]],coord[[k]]]-D[chris[[i,j,k]],coord[[ll]]]+Sum[chris[[i,k,s]]*chris[[s,ll,j]]-chris[[i,ll,s]]*chris[[s,k,j]],{s,1,n}],{i,1,n},{j,1,n},{k,1,n},{ll,1,n}];ricci=Table[Sum[riemann[[s,j,s,k]],{s,1,n}],{j,1,n},{k,1,n}];riccis=Sum[inversemetric[[i,j]]*ricci[[i,j]],{i,1,n},{j,1,n}];einstein=Table[ricci[[i,j]]-1/2*metric[[i,j]]*riccis,{i,1,n},{j,1,n}];


(* ::Input::Initialization:: *)
expricci=Simplify[Series[ricci,{po,0,1}]/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2iter,SphericalHarmonicY[l,m+3,th,ph]->ylm3iter},Assumptions->{l>m>0}];expricci>>"exp_e_ricci.txt";expeinstein=Simplify[Series[einstein,{po,0,1}]/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2iter,SphericalHarmonicY[l,m+3,th,ph]->ylm3iter},Assumptions->{l>m>0}];expeinstein>>"exp_e_eins.txt";


(* ::Input:: *)
(**)
(*(* Gauss-Bonnet energy-momentum tensor *)*)


(* ::Input::Initialization:: *)
riemanndown=Table[Sum[riemann[[s,j,k,ll]]*metric[[i,s]],{s,1,n}],{i,1,n},{j,1,n},{k,1,n},{ll,1,n}];gtensdown=Table[-riemanndown[[i,j,k,ll]]+metric[[i,k]]*ricci[[j,ll]]+metric[[j,ll]]*ricci[[i,k]]-metric[[i,ll]]*ricci[[j,k]]-metric[[j,k]]*ricci[[i,ll]]-1/2*(metric[[i,k]]*metric[[j,ll]]-metric[[j,k]]*metric[[i,ll]])*riccis,{i,1,n},{j,1,n},{k,1,n},{ll,1,n}];


(* ::Input::Initialization:: *)
dbps=Table[ddf*D[bph,coord[[i]]]*D[bph,coord[[j]]]+df*D[bph,coord[[i]],coord[[j]]]-df*Sum[chris[[i1,i,j]]*D[bph,coord[[i1]]],{i1,1,n}]/.{df->df+ddf*debph,ddf->ddf+dddf*debph},{i,1,n},{j,1,n}];tgbtensdown=Table[4*Sum[inversemetric[[i1,i2]]*inversemetric[[j1,j2]]*gtensdown[[i,i1,j,j1]]*dbps[[i2,j2]],{i1,1,n},{j1,1,n},{i2,1,n},{j2,1,n}],{i,1,n},{j,1,n}];tgbtenstr=Sum[inversemetric[[i1,i2]]*tgbtensdown[[i1,i2]],{i1,1,n},{i2,1,n}];sgbmn=Table[tgbtensdown[[i,j]]-1/2*metric[[i,j]]*tgbtenstr,{i,1,n},{j,1,n}];


(* ::Input::Initialization:: *)
exptgbmn=Simplify[Series[tgbtensdown,{po,0,1}]/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2iter,SphericalHarmonicY[l,m+3,th,ph]->ylm3iter},Assumptions->{l>m>0}];exptgbmn>>"exp_e_tgb.txt";expsgbmn=Simplify[Series[sgbmn,{po,0,1}]/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2iter,SphericalHarmonicY[l,m+3,th,ph]->ylm3iter},Assumptions->{l>m>0}];expsgbmn>>"exp_e_sgb.txt";


(* ::Input:: *)
(**)
(*(* scalar energy-momentum tensor *)*)


(* ::Input::Initialization:: *)
tphmn=Table[D[bph,coord[[i]]]*D[bph,coord[[j]]]-1/2*metric[[i,j]]*(Sum[inversemetric[[i1,j1]]*D[bph,coord[[i1]]]*D[bph,coord[[j1]]],{i1,1,n},{j1,1,n}]+bu)/.{bu->bu+dbu*debph},{i,1,n},{j,1,n}];tphtr=Sum[inversemetric[[i1,j1]]*tphmn[[i1,j1]],{i1,1,n},{j1,1,n}];sphmn=Table[tphmn[[i,j]]-1/2*metric[[i,j]]*tphtr,{i,1,n},{j,1,n}];


(* ::Input::Initialization:: *)
exptphmn=Simplify[Series[tphmn,{po,0,1}]/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2iter,SphericalHarmonicY[l,m+3,th,ph]->ylm3iter},Assumptions->{l>m>0}];exptphmn>>"exp_e_tph.txt";expsphmn=Simplify[Series[sphmn,{po,0,1}]/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2iter,SphericalHarmonicY[l,m+3,th,ph]->ylm3iter},Assumptions->{l>m>0}];expsphmn>>"exp_e_sph.txt";


(* ::Input:: *)
(**)
(*(* matter energy-momentum tensor *)*)


(* ::Input::Initialization:: *)
stensordown=Table[(ep+p)*u4veldown[[i]]*u4veldown[[j]]+1/2*(ep-p)*metric[[i,j]]/.{ep->epf[r]+deep,p->pf[r]+dep},{i,1,n},{j,1,n}];ttensordown=Table[(ep+p)*u4veldown[[i]]*u4veldown[[j]]+p*metric[[i,j]]/.{ep->epf[r]+deep,p->pf[r]+dep},{i,1,n},{j,1,n}];


(* ::Input::Initialization:: *)
exptmmn=Simplify[Series[ttensordown,{po,0,1}]/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2iter,SphericalHarmonicY[l,m+3,th,ph]->ylm3iter},Assumptions->{l>m>0}];exptmmn>>"exp_e_tm.txt";expsmmn=Simplify[Series[stensordown,{po,0,1}]/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2iter,SphericalHarmonicY[l,m+3,th,ph]->ylm3iter},Assumptions->{l>m>0}];expsmmn>>"exp_e_sm.txt";


(* ::Input::Initialization:: *)
ttensorup=Table[(ep+p)*u4velup[[i]]*u4velup[[j]]+p*inversemetric[[i,j]]/.{ep->epf[r]+deep,p->pf[r]+dep},{i,1,n},{j,1,n}];conseq=Table[Sum[D[ttensorup[[i1,i]],coord[[i1]]],{i1,1,n}]+Sum[chris[[j1,j1,j2]]*ttensorup[[j2,i]]+chris[[i,j1,j2]]*ttensorup[[j1,j2]],{j1,1,n},{j2,1,n}],{i,1,n}];


(* ::Input::Initialization:: *)
expconseq=Simplify[Series[conseq,{po,0,1}]/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2iter,SphericalHarmonicY[l,m+3,th,ph]->ylm3iter},Assumptions->{l>m>0}];expconseq>>"exp_e_coneq.txt";


(* ::Input:: *)
(**)
(*(* calculate the scalar equation *)*)


(* ::Input::Initialization:: *)
rgb=Sum[-inversemetric[[i1,i2]]*inversemetric[[i3,i4]]*inversemetric[[i5,i6]]*inversemetric[[i7,i8]]*gtensdown[[i1,i3,i5,i7]]*riemanndown[[i2,i4,i6,i8]],{i1,1,n},{i2,1,n},{i3,1,n},{i4,1,n},{i5,1,n},{i6,1,n},{i7,1,n},{i8,1,n}];scalareq=Sum[inversemetric[[i,j]]*D[D[bph,coord[[i]]],coord[[j]]],{i,1,n},{j,1,n}]-Sum[inversemetric[[i2,i3]]*chris[[i1,i2,i3]]*D[bph,coord[[i1]]],{i1,1,n},{i2,1,n},{i3,1,n}]+1/2*df*rgb-1/2*dbu/.{df->df+ddf*debph,dbu->dbu+ddbu*debph};


(* ::Input::Initialization:: *)
expscalareq=Simplify[Series[scalareq,{po,0,1}]/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2iter,SphericalHarmonicY[l,m+3,th,ph]->ylm3iter},Assumptions->{l>m>0}];expscalareq>>"exp_e_scalareq.txt";


(* ::Input:: *)
(**)
(**)
(*(* II. remove the angular dependence to get equations in variables of t and r *)*)


(* ::Input:: *)
(*(* defind the convariant basis *)*)


(* ::Input::Initialization:: *)
basisv={D[SphericalHarmonicY[l,m,th,ph],th],D[SphericalHarmonicY[l,m,th,ph],ph]};basist={D[D[SphericalHarmonicY[l,m,th,ph],th],th],D[D[SphericalHarmonicY[l,m,th,ph],ph],th]-Cot[th]*D[SphericalHarmonicY[l,m,th,ph],ph],D[D[SphericalHarmonicY[l,m,th,ph],ph],ph]+Sin[th]*Cos[th]*D[SphericalHarmonicY[l,m,th,ph],th]};basisv=Simplify[basisv/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2iter},Assumptions->{l>m>0}];basist=Simplify[basist/.{SphericalHarmonicY[l,m,th,ph]->ylm,SphericalHarmonicY[l,m+1,th,ph]->ylm1,SphericalHarmonicY[l,m+2,th,ph]->ylm2iter},Assumptions->{l>m>0}];


(* ::Input:: *)
(*(* read the data from part I *)*)


(* ::Input::Initialization:: *)
ricci=<<"exp_e_ricci.txt";sgbmn=<<"exp_e_sgb.txt";sphmn=<<"exp_e_sph.txt";smmn=<<"exp_e_sm.txt";scalareq=<<"exp_e_scalareq.txt";
conseq=<<"exp_e_coneq.txt";
einstein=<<"exp_e_eins.txt";tgbmn=<<"exp_e_tgb.txt";tphmn=<<"exp_e_tph.txt";tmmn=<<"exp_e_tm.txt";ricci0=Coefficient[ricci,po,0];riccip=Coefficient[ricci,po,1];sgbmn0=Coefficient[sgbmn,po,0];sgbmnp=Coefficient[sgbmn,po,1];sphmn0=Coefficient[sphmn,po,0];sphmnp=Coefficient[sphmn,po,1];smmn0=Coefficient[smmn,po,0];smmnp=Coefficient[smmn,po,1];scalareq0=Coefficient[scalareq,po,0];scalareqp=Coefficient[scalareq,po,1];
conseq0=Coefficient[conseq,po,0];conseqp=Coefficient[conseq,po,1];einstein0=Coefficient[einstein,po,0];einsteinp=Coefficient[einstein,po,1];tgbmn0=Coefficient[tgbmn,po,0];tgbmnp=Coefficient[tgbmn,po,1];tphmn0=Coefficient[tphmn,po,0];tphmnp=Coefficient[tphmn,po,1];tmmn0=Coefficient[tmmn,po,0];tmmnp=Coefficient[tmmn,po,1];


(* ::Input:: *)
(*(* 0th-order equations *)*)


(* ::Input::Initialization:: *)
eeq0=Simplify[Table[einstein0[[i,j]]-tgbmn0[[i,j]]-tphmn0[[i,j]]-ka*tmmn0[[i,j]],{i,1,n},{j,1,n}]];eq0={eeq0[[1]],eeq0[[2]],eeq0[[3]],eeq0[[4]],conseq0[[1]],scalareq0};eq0>>"eqs_order0.txt";


(* ::Input:: *)
(*listeq0:=Table[If[UnsameQ[eeq0[[i,j]],0],{ToString[Eq[i,j]],eeq0[[i,j]]}],{i,1,n},{j,1,n}];TableForm[Partition[DeleteCases[Flatten[listeq0],Null],2],TableSpacing->{2,2}]*)


(* ::Input:: *)
(*(* 1st-order equations *)*)


(* ::Input::Initialization:: *)
eeq1=Table[einsteinp[[i,j]]-tgbmnp[[i,j]]-tphmnp[[i,j]]-ka*tmmnp[[i,j]],{i,1,n},{j,1,n}];deeqe0=Coefficient[eeq1[[1,1]],ylm];
deeqe1=Coefficient[eeq1[[1,2]],ylm];deeqe2=Coefficient[eeq1[[1,4]],ylm]/(I*m);deeqe3=Coefficient[eeq1[[2,2]],ylm];deeqe4=Coefficient[eeq1[[2,4]],ylm]/(I*m);deeqe6=FullSimplify[Coefficient[eeq1[[3,4]],ylm1],Assumptions->{l>m>0}]/( I E^(-I ph) m \[Sqrt]((l-m)(l+m+1)));test=Simplify[eeq1[[3,3]]-deeqe6*basist[[1]],Assumptions->{l>m>0}];deeqe5=Coefficient[test,ylm]/r^2;descalareq=Coefficient[scalareqp,ylm];deconseq1=Coefficient[conseqp[[1]],ylm];deconseq2=Coefficient[conseqp[[2]],ylm];deconseq3=(r^2*(Sin[th])^2*Coefficient[conseqp[[4]],ylm])/(I*m);eq1={deeqe0,deeqe1,deeqe2,deeqe3,deeqe4,deeqe5,deeqe6,deconseq1,deconseq2,deconseq3,descalareq};eq1>>"eqs_even_order1.txt";
