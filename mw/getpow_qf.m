function pow = getpow_qf(om0,om1,e)
om0i=inv(om0);
om1i=inv(om1);
ch_om0 = chol(om0);
ch_om1 = chol(om1);
ch_om0i = chol(om0i);
ch_om1i = chol(om1i);
ho=ch_om1i*ch_om0';
ha=ch_om0i*ch_om1';
qe = sum(e.^2,1)';
ya_o = ho*e;
yo_a = ha*e;
qa_o = sum(ya_o.^2,1)';
qo_a = sum(yo_a.^2,1)';
lr_o = qe./qa_o;
lr_a = qo_a./qe;
cv = prctile(lr_o,95);
pow = mean(lr_a>cv);

end