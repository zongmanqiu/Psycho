* Encoding: UTF-8.
PRESERVE.
set printback=off.
set LENGTH=256.
set mxloops=1000000000.
/** inverse cdf of normal distribution **/.
define @ppnd16(idfp=!charend('/'))
compute invnp=!idfp.
compute invna0=3.3871328727963666080.
compute invna1=133.14166789178437745.
compute invna2=1971.5909503065514427.
compute invna3=13731.693765509461125.
compute invna4=45921.953931549871457.
compute invna5=67265.770927008700853.
compute invna6=33430.575583588128105.
compute invna7=2509.0809287301226727.
compute invnb1=42.313330701600911252.
compute invnb2=687.18700749205790830.
compute invnb3=5394.1960214247511077.
compute invnb4=21213.794301586595867.
compute invnb5=39307.895800092710610.
compute invnb6=28729.085735721942674.
compute invnb7=5226.4952788528545610.
compute invnc0=1.42343711074968357734.
compute invnc1=4.63033784615654529590.
compute invnc2=5.76949722146069140550.
compute invnc3=3.64784832476326738258.
compute invnc4=1.27045825245236838258.
compute invnc5=0.241780725177450611770.
compute invnc6=0.0227238449892691845833.
compute invnc7=0.000774545014278341407640.
compute invnd1=2.05319162663775882187.
compute invnd2=1.67638483018380384940.
compute invnd3=0.689767334985100004550.
compute invnd4=0.148103976427480074590.
compute invnd5=0.0151986665636164571966.
compute invnd6=0.000547593808499534494600.
compute invnd7=0.105075007164441684324.
compute invne0=6.65790464350110377720.
compute invne1=5.46378491116411436990.
compute invne2=1.78482653991729133580.
compute invne3=0.296560571828504891230.
compute invne4=0.0265321895265761230930.
compute invne5=0.00124266094738807843860.
compute invne6=0.0000271155556874348757815.
compute invne7=0.000000201033439929228813265.
compute invnf1=0.599832206555887937690.
compute invnf2=0.136929880922735805310.
compute invnf3=0.0148753612908506148525.
compute invnf4=0.000786869131145613259100.
compute invnf5=0.0000184631831751005468180.
compute invnf6=0.000000142151175831644588870.
compute invnf7=0.204426310338993978564.
compute invnq=invnp-0.5.
do if (abs(invnq) <= 0.425).
    compute invnr=0.180625-invnq**2.
    compute ppnd16fz=invnq*(((((((invna7*invnr+invna6)*invnr+invna5)*invnr+
        invna4)*invnr+invna3)*invnr+invna2)*invnr+invna1)*invnr+invna0).
    compute ppnd16fm=(((((((invnb7*invnr+invnb6)*invnr+invnb5)*invnr+invnb4)*
        invnr+invnb3)*invnr+invnb2)*invnr+invnb1)*invnr+1).
    compute ppnd16=ppnd16fz/ppnd16fm.
    else.
    do if (invnq<=0).
        compute invnr=invnp.
        ELSE.  
        compute invnr=1-invnp.
    end if.
    do if (invnr<=0).
        compute ppnd16=0.
        else.
        compute invnr=sqrt(-ln(invnr)).
        do if (invnr<=5).
            compute invnr=invnr-1.6.
            compute ppnd16fz=(((((((invnc7*invnr+invnc6)*invnr+invnc5)*invnr+
                invnc4)*invnr+invnc3)*invnr+invnc2)*invnr+invnc1)*invnr+invnc0).
            compute ppnd16fm=(((((((invnd7*invnr*0.00000001+invnd6)*invnr+
                invnd5)*invnr+invnd4)*invnr+invnd3)*invnr+invnd2)*invnr+invnd1)*invnr+1).
            compute ppnd16=ppnd16fz/ppnd16fm.
            else.
            compute invnr=invnr-5. 
            compute ppnd16fz=(((((((invne7*invnr+invne6)*invnr+invne5)*invnr+
                invne4)*invnr+invne3)*invnr+invne2)*invnr+invne1)*invnr+invne0).
            compute ppnd16fm=(((((((invnf7*invnr*0.00000000000001+invnf6)*
                invnr+invnf5)*invnr+invnf4)*invnr+invnf3)*invnr+invnf2)*invnr+invnf1)*invnr+1).
            compute ppnd16=ppnd16fz/ppnd16fm.
        end if.
        do if (invnq<0).
            compute ppnd16=-ppnd16.
        end if.
    end if.
end if.
!enddefine.
/** Bias-corrected CI **/.
define @bcci(bcciy=!charend('/')/bccih=!default(97.5) !charend('/')/bcciraw=!charend('/'))
compute bccixx=!bcciy.
compute bccip3=!bccih.
@ppnd16 idfp=bccip3.
compute bccirz=ppnd16.
compute bccin=nrow(bccixx).
compute bcciul=99.94-98/bccin.
compute bccill=0.01+99/bccin.
compute bccinc=0.
loop #bc=1 to bccin.
    do if (bccixx(#bc,1)<!bcciraw).
        compute bccinc=bccinc+1.
    end if.
end loop.
compute bccinp=bccinc/(bccin+1).
@ppnd16 idfp=bccinp.
compute bccimz=ppnd16.
compute bccilz=2*bccimz-bccirz.
compute bccihz=2*bccimz+bccirz.
compute bccilp=cdfnorm(bccilz)*100.
compute bccihp=cdfnorm(bccihz)*100.
do if (bccihp>bcciul).
    compute bccihp=bcciul.
    do if (csum(ero)<>12).
        compute ero={ero;12}.
    end if.
end if.
do if (bccihp<bccill).
    compute bccihp=bccill.
    do if (csum(ero)<>12).
        compute ero={ero;12}.
    end if.
end if.
do if (bccilp<bccill).
    compute bccilp=bccill.
    do if (csum(ero)<>12).
        compute ero={ero;12}.
    end if.
end if.
do if (bccilp>bcciul).
    compute bccilp=bcciul.
    do if (csum(ero)<>12).
        compute ero={ero;12}.
    end if.
end if.
@prcn prcny=bccixx/prcnp=bccilp.
compute bccilp=percnv.
@prcn prcny=bccixx/prcnp=bccihp.
compute bccilh=percnv.
!enddefine.
/** Percentile **/.
define @prcn(prcny=!charend('/')/prcnp=!charend('/'))
compute percnp=!prcnp.
compute percny=!prcny.
compute percnn=nrow(percny).
compute percnnn=(percnn+1)*percnp/100.
compute percni=trunc(percnnn).
compute percni2=percni+1.
compute percnd=percnnn-percni.
do if (percnd =0).
    compute  percnv=percny(percni,1).
    else.
    compute  percnv=percnd*percny(percni2,1)+(1-percnd)*percny(percni,1).
end if.
!enddefine.
/** Bootstrap's output **/.
define @btsot(btrawd=!charend('/')/btalld=!charend('/')/btcitp=!default(1) !charend('/')
    /btnewnm=!charend('/')/btcih=!default(97.5) !charend('/')/btcil=!default(2.5) !charend('/'))
    compute btcinvar=ncol(!btrawd).
    compute btcin=nrow(!btalld).
    compute !btnewnm=make(6,btcinvar,-999).
    compute !btnewnm(1,:)=!btrawd.
    compute !btnewnm(2,:)=csum(!btalld)/btcin.
    loop #btc=1 to btcinvar.
        compute btmpdt=!btalld(:,#btc).
        compute btmpdt2=make(btcin,1,!btnewnm(2,#btc)).
        compute !btnewnm(4,#btc)=(csum((btmpdt-btmpdt2)&**2)/(btcin-1))**0.5.
        @prcn prcny=btmpdt/prcnp=50.
        compute !btnewnm(3,#btc)=percnv.
        do if (!btcitp=1).
            @prcn prcny=btmpdt/prcnp=!btcil.
            compute !btnewnm(5,#btc)=percnv. 
            @prcn prcny=btmpdt/prcnp=!btcih.
            compute !btnewnm(6,#btc)=percnv.
            else.    
            @bcci bcciy=btmpdt/bccih=bcrcih/bcciraw=!btrawd(1,#btc).
            compute !btnewnm(5,#btc)=bccilp. 
            compute !btnewnm(6,#btc)=bccilh. 
        end if.
    end loop.
!enddefine.
/** Data sorting **/.
define @sortv(sortv=!charend('/')/sorttp=!default(1) !charend('/'))
    compute sortdt=!sortv.
    compute sorttp=!sorttp.    
    compute sort2=grade(sorttp&*sortdt).
    compute sortdt(sort2)=sortdt.
    compute !sortv=sortdt.
    release sortdt.
!enddefine.
/** Check constant and correlation **/.
define @candr(crdata=!charend('/'))
    compute crdata=!crdata.
    compute crnvar=ncol(crdata).
    compute candrn=0.
    loop #i=1 to crnvar.
        compute tdata=crdata(1:ncase,#i).
        compute tmax=cmax(tdata).
        compute tmin=cmin(tdata).
        do if (tmax=tmin).
            compute candrn=1.
        end if.
    end loop.
    do if (candrn=0).
        compute tpar1=1/(ncase-1).
        compute corvcv=tpar1*(sscp(crdata)-((t(csum(crdata))*csum(crdata))/ncase)).
        compute cord=inv(mdiag(sqrt(diag(corvcv)))).
        compute corcr=cord*corvcv*cord.
        loop #ckcr1=1 to crnvar.
            loop #ckcr2 =1 to crnvar.
                do if (#ckcr1 <> #ckcr2).
                    do if (corcr(#ckcr1,#ckcr2)=1).
                        compute candrn=2.
                    end if.
                end if.
            end loop.
        end loop.
    end if.
    do if (candrn=0).
        compute ckeval=eval(corcr).
        do if (cmin(ckeval)<0.00001).
            compute candrn=3.
        end if.
    end if.
    release crdata,crnvar,tdata,tmax,tmin,tpar1,corvcv,cord,corcr,ckeval.
!enddefine.
/**htmt analysis.
define @htmt(htmtdt=!charend('/'))
compute varvcv=(sscp(!htmtdt)-(sscp(csum(!htmtdt))/ncase))/(ncase-1)./*var-covariance matrix.
compute varsd=sqrt(diag(varvcv)).  /*std.dev.
compute varcorr=abs(inv(mdiag(varsd))*varvcv*inv(mdiag(varsd))).  /*absolute correlations. 
compute htmtvcv=make(nf,nf,-999).  /*htmt matrix.
loop #i1=1 to nf.
    loop #i2=1 to nf.
        do if (#i1 = #i2). /*diagonal elements.
            compute rstar=rsum(fnum(1,1:#i1))+1.
            compute rend=rsum(fnum(1,1:(#i1+1))).
            compute cstar=rstar.
            compute cend=rend.
            compute htmtvcv(#i1,#i1)=(csum(rsum(varcorr(rstar:rend,cstar:cend)))-
                fnum(1,(#i1+1)))/(fnum(1,(#i1+1))**2-fnum(1,(#i1+1))).
            else if (#i1<#i2). /*non-diagonal elements.
            compute rstar=rsum(fnum(1,1:#i1))+1.
            compute rend=rsum(fnum(1,1:(#i1+1))).
            compute cstar=rsum(fnum(1,1:#i2))+1.
            compute cend=rsum(fnum(1,1:(#i2+1))).
            compute htmtvcv(#i1,#i2)=csum(rsum(varcorr(rstar:rend,cstar:cend)))
                /fnum(1,(#i1+1))/fnum(1,(#i2+1)).
            compute htmtvcv(#i2,#i1)=htmtvcv(#i1,#i2). 
        end if.
    end loop.
end loop.
compute htmtsd=sqrt(diag(htmtvcv)).
compute htmtcorr=inv(mdiag(htmtsd))*htmtvcv*inv(mdiag(htmtsd)). /*htmt matrix.
!enddefine. 
define Psycho (var= !charend('/')
/fnum=!charend('/')
/dcms=!default(f12.8) !charend('/')
/bootn=!default(1000) !charend('/')
/bootcim=!default(1) !charend('/')
/cinum=!default(95) !charend('/'))
matrix.
    get raw//file=*/var=!var /name=nmall/missing=omit.
    compute fnum={0,!fnum}. /*n of items in each factor. 
    compute sumf=rsum(fnum). /*total n of items.
    compute nf=ncol(fnum)-1. /*n of factors.
    compute ncase=nrow(raw). /*n of cases.
    compute nvar=ncol(raw). /*n of variables.
    compute bootn=!bootn. /*n of bootstrap.
    !let !tt=('*********************************************'+
        '********************************************')
    !let !tt=!quote(!tt)
    print /title=!tt/format=a8.
    print /title="Psycho version 1.0(Beta)"/format=a8/space=0.
    print /title='Copyright (C) 2022 by Qiu Zongman'/format=a8.
    print /title=!tt/format=a8.
    compute fnm={"F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12","F13","F14",
        "F15","F16","F17","F18","F19","F20","F21","F22","F23","F24","F25","F26","F27","F28",
        "F29","F30","F31","F32","F33","F34","F35","F36","F37","F38","F39","F40","F41","F42",
        "F43","F44","F45","F46","F47","F48","F49","F50"}.
    print {ncase,nvar,nf,!bootn,!cinum}/title='Number of ...'
        /cname={'cases','items','factors','Boot_n','CI'}/format=!dcms.
    do if (!bootcim=1).
        print /title ='Bootstrap CI: Percentiles'.
        else.
        print /title ='Bootstrap CI: Bias-Corrected'.
    end if.
    do if (sumf=nvar). 
        /** factors by items.
        compute fimax=rmax(fnum).
        compute fnmp=make(nf,fimax," ").
        compute floop=0.
        loop i0=1 to nf.
            compute tmp1=fnum(1,(i0+1)).
            loop #i02=1 to tmp1.
                compute floop=floop+1.
                compute fnmp(i0,#i02)=nmall(1,floop).
            end loop.
        end loop.
        compute fnmp2={fnm(1,1:nf),"       ."}.
        compute fnmp={t(fnm(1,1:nf)),make(nf,1,"by"),fnmp}.
        print fnmp/title='>>>>> Factors by items'/format=a8.
        @htmt htmtdt=raw.
        print htmtcorr//title='>>>>> HTMT matrix '/cname=fnmp2/rname=fnmp2/format=!dcms.
        compute npair=nf*(nf-1)/2.
        compute htmtpair=make(npair,(bootn+3),-999).
        compute floop=0.
        loop #i1=1 to nf. /*matrix to column.
            loop #i2=1 to nf.
                do if (#i1<#i2).
                    compute floop=floop+1.
                    compute htmtpair(floop,1)=#i1.
                    compute htmtpair(floop,2)=#i2.
                    compute htmtpair(floop,3)=htmtcorr(#i1,#i2).
                end if.
            end loop.
        end loop.
        /*bootstrap .
        do if (!bootn>0).
            compute btsmpl=raw.
            compute btero=0.
            compute btsmpln=0.
            loop #bt1=1 to !bootn.
                loop #bt2=1 to ncase.
                    compute btsmpl(#bt2,:)=raw((trunc(uniform(1,1)*ncase)+1),:).
                end loop.
                @candr crdata=btsmpl.
                do if (candrn>0).
                    compute btero=btero+1.
                    else.
                    @htmt htmtdt=btsmpl.
                    compute btsmpln=btsmpln+1.
                    compute floop=0.
                    loop #i3=1 to nf.
                        loop #i4=1 to nf.
                            do if (#i3<#i4).
                                compute floop=floop+1.
                                compute htmtpair(floop,(3+btsmpln))=htmtcorr(#i3,#i4).
                            end if.
                        end loop.
                    end loop.
                end if.
            end loop.
            do if (btero>0).
                compute btero2=bootn-btero.
                loop #i1=1 to btero.
                    compute htmtpair(:,(btsmpln+4-#i1))=
                    htmtpair(:,(trunc(uniform(1,1)*(bootn-btero))+4))).
                end loop.
            end if .
            compute htmt1st=htmtpair(:,1:2).
            compute htmt2nd=htmtpair(:,4:(btsmpln+3)).
            compute htmt2nd=t(htmt2nd).
            loop #i5=1 to npair.
                @sortv sortv=htmt2nd(:,#i5).
            end loop.
            compute cih=50+!cinum/2.
            compute cil=100-cih.
            @btsot btrawd=t(htmtpair(:,3))/btalld=htmt2nd/btcitp=!bootcim/btnewnm=btsmpl
                /btcih=cih/btcil=cil.
            compute htmtpair={htmt1st,t(btsmpl)}.
            print htmtpair/title='>>>>> Bootstrap'/cname={'Factor A','Factor B','HTMT',
                'BootMean','BootMdn','BootSE','BootLLCI','BootULCI'}/rname={"Results:"}
                /format=!dcms.
            do if (btero>0).
                print /title='************************************ Errors and Notes '+
                '***********************************'/format=a8.
                print btero/title'Note   : Due to estimation problems, some bootstrap '
                +'samples had to be replaced.' /rname={'Size   :'}/format=!dcms.
            end if.
        end if.
    else.
        print /title='************************************ Errors and Notes '+
        '***********************************'/format=a8.
        print /title='Error  : Invalid factors.'/format=a8.
    end if.
end matrix.
!enddefine.
Psycho var=%%var%%
/fnum=%%fnum%%
/dcms=%%dcms%%
/bootn=%%bootn%%
/bootcim=%%bootcim%%
/cinum=%%cinum%%.
RESTORE.
