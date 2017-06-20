unit drag;
interface
Uses Dialogs,Math,SysUtils;
const korn=200;
type
  nla=array[1..100,1..100] of real;
  nly=array[1..100] of real;
  nlaux=array[1..100,1..100] of real;

var
  aprl,aprr             : array[1..100,1..200] of real;
  a,aap                 : nla;
  y,yy,y0,prmt,eql,dery : nly;
  ytrue,R               : nly;
  aux                   : nlaux;
  ii,i,n,ihlf,j0,jk,jc  : integer;
  ntrue,ng,nb           : integer;
  d,direk               : real;
  Kn                    : real;
  gg                    : boolean;
  mhu,Zu,Zf             : array [1..100] of integer;


var sng_real,sng_imag:nly;
var ys_real,ys_imag:nla;
var yr_baz,yi_baz:Real;

      procedure f1(var {eql,}y:nly);
      procedure f2(var aa:nla;var y:nly);
      procedure nldet(n:integer;var a:nla;var d:real);
      procedure nleqc(ndim:integer;var x,direk:real;var y,yy,y0,{eql,}dery,prmt:nly;var aa,aap:nla);
      procedure nleqd(ndim,ihlf:integer;var x,direk:real;var y,yy,{eql,}prmt:nly;var aa,aap:nla);
      procedure nlmain(ndim:integer;var ihlf:integer;var y0,dery,eql,y,yy,prmt:nly;var aa,aap:nla;var aux:nlaux;var direk:real);
      procedure dragmain;
      function fi(i:integer;y1:real;y2:real;f1:real;f2:real):real;
implementation

//{$R *.DFM}

procedure f1(var {eql,}y:nly);
var kJac,k,j,l:Integer;
var tmp1,tmp2:double;
begin

 {This procedure calculates values of the equations for y[i]  }
 {User fills this procedure before the compilation of the unit}
//���� ��� P
j:=1;
for k:=1 to ntrue do
 if((mhu[k]<>0)and(mhu[k]<>4))then
   begin
   ytrue[k]:=y[j];
   j:=j+1;
   end;
for k:=1 to ntrue do
 if(mhu[k]<>0)then
   begin
   ytrue[k+ntrue]:=y[j];
   j:=j+1;
   end;
for k:=1 to (n div 2)-1 do
 R[k]:=y[k+(n div 2)];

//R[(n div 2)]:=-1;
Kn:=y[n];

kJac := 1;
for k := 1 to ntrue do
 if (mhu[k] <> 0) then
 begin
//  Nebal[kJac]=0;
  tmp1 := 0;
  tmp2 := 0;
  tmp1 := -sng_real[k];
  tmp2 :=tmp2+ ys_real[k,k]*ytrue[k]*ytrue[k];
  for j := 1 to ntrue do
    if (j <>k) then
    begin
     tmp2 :=tmp2+ ytrue[k]*ytrue[j]*(ys_real[k,j]*(cos(ytrue[ntrue+k])*cos(ytrue[ntrue+j])+sin(ytrue[ntrue+k])*sin(ytrue[ntrue+j]))-ys_imag[k,j]*(cos(ytrue[ntrue+k])*sin(ytrue[ntrue+j])-sin(ytrue[ntrue+k])*cos(ytrue[ntrue+j])));
    end;
//  aa[kJac,n-1] := tmp2;
  eql[kJac] := tmp1+tmp2*Kn;
 // Nebal[kJac]=tmp1+tmp2;
  kJac:=kJac+1;
 end;

//���� ��� Q
kJac := 1;
for k := 1 to ntrue do
if ((mhu[k] <> 0) and (mhu[k] <> 4)) then
 begin
//  Nebal[kJac+n-nb]=0;
  tmp2 := 0;
  tmp1 := 0;
  tmp1 := -sng_imag[k];
  tmp2 :=tmp2+ ys_imag[k,k]*ytrue[k]*ytrue[k];
  for j := 1 to ntrue do
    if (j <> k) then
    begin
    tmp2 :=tmp2+ ytrue[k]*ytrue[j]*(ys_imag[k,j]*(cos(ytrue[ntrue+k])*cos(ytrue[ntrue+j])+sin(ytrue[ntrue+k])*sin(ytrue[ntrue+j]))+ys_real[k,j]*(cos(ytrue[ntrue+k])*sin(ytrue[ntrue+j])-sin(ytrue[ntrue+k])*cos(ytrue[ntrue+j])));
    end;
//  aa[kJac+n-nb,n-1] := tmp2;
  eql[kJac+ntrue-nb] := tmp1+tmp2*Kn;
//  Nebal[kJac+n-nb]=tmp1+tmp2;
  kJac:=kJac+1;
 end;

//���� ��� Vp
kJac := 1;
for k := 1 to ntrue do
 if (mhu[k] <> 0) then
 begin
  tmp1 := 0;
  tmp2 := 0;
  for l := 1 to ntrue do
   if((mhu[l] <> 0) and (mhu[l] <> 4)) then
   begin
    if (k=l) then
     begin
      tmp1 :=tmp1+ 2.*ys_real[k,k]*ytrue[k]*R[Zu[k]];
      for j := 1 to ntrue do
       if (j <> k) then
         tmp1 :=tmp1+ R[Zu[l]]*ytrue[j]*fi(1,ys_real[k,j],-ys_imag[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
     end
    else
      tmp1 :=tmp1+ R[Zu[l]]*ytrue[k]*fi(1,ys_real[k,l],-ys_imag[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
   end;
  for l := 1 to ntrue do
   if (mhu[l] <> 0) then
   begin
    if (k=l) then
     begin
      for j := 1 to ntrue do
       if (j <> k) then
         tmp2 :=tmp2+ R[Zf[l]]*ytrue[k]*ytrue[j]*fi(2,ys_real[k,j],-ys_imag[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
     end
    else
      tmp2 :=tmp2+ R[Zf[l]]*ytrue[k]*ytrue[l]*fi(3,ys_real[k,l],-ys_imag[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
   end;
// aa[(n div 2)+kJac,n] := tmp1+tmp2;
 eql[(n div 2)+kJac] := (tmp1*Kn+tmp2*Kn);
 kJac:=kJac+1;
 end;

//Vq
kJac := 1;
for k := 1 to ntrue do
 if ((mhu[k] <> 0)and (mhu[k]<>4)) then
 begin
  tmp1 := 0;
  tmp2 := 0;
  for l := 1 to ntrue do
   if((mhu[l] <> 0) and (mhu[l] <> 4)) then
   begin
    if (k=l) then
     begin
      tmp1 :=tmp1+ 2.*ys_imag[k,k]*ytrue[k]*R[Zu[k]];
      for j := 1 to ntrue do
       if (j <> k) then
         tmp1 :=tmp1+ R[Zu[l]]*ytrue[j]*fi(1,ys_imag[k,j],ys_real[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
     end
    else
      tmp1 :=tmp1+ R[Zu[l]]*ytrue[k]*fi(1,ys_imag[k,l],ys_real[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
   end;
  for l := 1 to ntrue do
   if (mhu[l] <> 0) then
   begin
    if (k=l) then
     begin
      for j := 1 to ntrue do
       if (j <> k) then
         tmp2 :=tmp2+ R[Zf[l]]*ytrue[k]*ytrue[j]*fi(2,ys_imag[k,j],ys_real[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
     end
    else
      tmp2 :=tmp2+ R[Zf[l]]*ytrue[k]*ytrue[l]*fi(3,ys_imag[k,l],ys_real[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
   end;
// aa[(n div 2)+kJac+ntrue-nb,n] := tmp1+tmp2;
 eql[(n div 2)+kJac+ntrue-nb] := (tmp1*Kn+tmp2*Kn);
 kJac:=kJac+1;
 end;
{

kJac:=1;
for k:=1 to ntrue do
begin
 if (mhu[k]<>0) then
  begin
   tmp1:=0;
   tmp1:=-sng_real[k]+ys_real[k,k]*ytrue[k]*ytrue[k];
   for j:=1 to ntrue do
    begin
     if (j<>k) then
      begin
       tmp1:=tmp1+ytrue[k]*ytrue[j]*(ys_real[k,j]*(cos(ytrue[k+ntrue])*cos(ytrue[j+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[j+ntrue]))-ys_imag[k,j]*(cos(ytrue[k+ntrue])*sin(ytrue[j+ntrue])-sin(ytrue[k+ntrue])*cos(ytrue[j+ntrue])));
      end;
    end;
  eql[kJac]:=tmp1;
  kJac:=kJac+1;
  end;
end;

//���� ��� Q
kJac:=1;
for k:=1 to ntrue do
begin
 if ((mhu[k]<>0)and(mhu[k]<>4)) then
  begin
   tmp2:=0;
   tmp2:=-sng_imag[k]+ys_imag[k,k]*ytrue[k]*ytrue[k];
   for j:=1 to ntrue do
    begin
     if (j<>k) then
      begin
       tmp2:=tmp2+ytrue[k]*ytrue[j]*(ys_imag[k,j]*(cos(ytrue[k+ntrue])*cos(ytrue[j+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[j+ntrue]))+ys_real[k,j]*(cos(ytrue[k+ntrue])*sin(ytrue[j+ntrue])-sin(ytrue[k+ntrue])*cos(ytrue[j+ntrue])));
      end;
     end;
  eql[kJac+ntrue-nb]:=tmp2;
  kJac:=kJac+1;
  end;
end;
}

end;

procedure f2(var aa:nla;var y:nly);
var
  dwPdr,dwQdr,dwPdi,dwQdi,tmp1,tmp2:Real;
var i,k,j,l,kJac,iJac,inttemp1,inttemp2:Integer;
//tmp_i,tmp_r:Real;

begin
j:=1;
for k:=1 to ntrue do
 if((mhu[k]<>0)and(mhu[k]<>4))then
   begin
   ytrue[k]:=y[j];
   j:=j+1;
   end;
for k:=1 to ntrue do
 if(mhu[k]<>0)then
   begin
   ytrue[k+ntrue]:=y[j];
   j:=j+1;
   end;
for k:=1 to (n div 2)-1 do
 R[k]:=y[k+(n div 2)];

//R[(n div 2)]:=-1;
Kn:=y[n];
//���� ��� dP/dU

for k:=1 to n do
 for j:=1 to n do
   aa[k,j]:=0;

kJac:=1;
for k := 1 to ntrue do
 if (mhu[k] <> 0) then
 begin
 iJac := 1;
 for i := 1 to ntrue do
 if ((mhu[i] <> 0) and (mhu[i] <> 4)) then
 begin
  dwPdr := 0;
  if (i=k) then
  begin
  dwPdr :=dwPdr+ 2.*ys_real[k,k]*ytrue[k];
  for j := 1 to ntrue do
    begin
        if (j <> k) then
         begin
          dwPdr:=dwPdr+ ytrue[j]*(ys_real[k,j]*(cos(ytrue[k+ntrue])*cos(ytrue[j+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[j+ntrue]))-ys_imag[k,j]*(cos(ytrue[k+ntrue])*sin(ytrue[j+ntrue])-sin(ytrue[k+ntrue])*cos(ytrue[j+ntrue])));
         end;
    end;
  end
  else
  begin
   dwPdr:=dwPdr+ ytrue[k]*(ys_real[k,i]*(cos(ytrue[k+ntrue])*cos(ytrue[i+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[i+ntrue]))-ys_imag[k,i]*(cos(ytrue[k+ntrue])*sin(ytrue[i+ntrue])-sin(ytrue[k+ntrue])*cos(ytrue[i+ntrue])));
  end;
 dwPdr := dwPdr*Kn;
 aa[kJac,iJac] := dwPdr;
 iJac:=iJac+1;
 end;
 kJac:=kJac+1;
 end;

iJac:=1;
kJac:=1;

//���� ��� dQ/dU
for k :=1 to ntrue do
 if ((mhu[k] <> 0) and (mhu[k] <> 4)) then
 begin
 iJac:=1;
 for i :=1 to ntrue do
 if ((mhu[i] <> 0) and (mhu[i] <> 4)) then
 begin
  dwQdr :=0;
  if (i=k) then
  begin
  dwQdr :=dwQdr+2.*ys_imag[k,k]*ytrue[k];
  for j :=1 to ntrue do
    begin
        if (j <> k) then
         begin
          dwQdr:= dwQdr+ytrue[j]*(ys_imag[k,j]*(cos(ytrue[k+ntrue])*cos(ytrue[j+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[j+ntrue]))+ys_real[k,j]*(cos(ytrue[k+ntrue])*sin(ytrue[j+ntrue])-sin(ytrue[k+ntrue])*cos(ytrue[j+ntrue])));
         end;
    end;
  end
  else
  begin
   dwQdr:= dwQdr+ytrue[k]*(ys_imag[k,i]*(cos(ytrue[k+ntrue])*cos(ytrue[i+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[i+ntrue]))+ys_real[k,i]*(cos(ytrue[k+ntrue])*sin(ytrue[i+ntrue])-sin(ytrue[k+ntrue])*cos(ytrue[i+ntrue])));
  end;
 dwQdr :=dwQdr* Kn;
 aa[kJac+ntrue-nb,iJac] := dwQdr;
 iJac:=iJac+1;
 end;
 kJac:=kJac+1;
 end;

iJac:=1;
kJac:=1;

//***********************
//���� ��� dVp/dU

for k:=1 to ntrue do
 if (mhu[k] <> 0) then
 begin
 iJac:=1;
 for i:=1 to ntrue do
 if ((mhu[i] <> 0) and (mhu[i] <> 4)) then
 begin
 dwPdr := 0;
 dwQdr := 0;
 if (i=k) then
  begin
   if (mhu[k] <> 4) then
   begin
    dwPdr:=dwPdr+ 2.*ys_real[k,k]*R[(Zu[k])];
    for j:=1 to ntrue do
     if (j <> k) then
      begin
      dwQdr:=dwQdr+ R[(Zf[k])]*ytrue[j]*fi(2,ys_real[k,j],-ys_imag[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
      end;
    for l:=1 to ntrue do
     if ((l <> k) and (mhu[l] <> 0)) then
      begin
       dwQdr:=dwQdr+ R[Zf[l]]*ytrue[l]*fi(3,ys_real[k,l],-ys_imag[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
       if (mhu[l] <> 4) then
        dwPdr:=dwPdr+ R[Zu[l]]*fi(1,ys_real[k,l],-ys_imag[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
      end;
   end
   else //��� ������������� ����
   begin
    for j:=1 to ntrue do
     if (j <> k) then
      dwQdr:=dwQdr+ R[Zf[k]]*ytrue[j]*fi(2,ys_real[k,j],-ys_imag[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
    for l:=1 to ntrue do
    begin
     if ((mhu[l] <> 0) and (mhu[l] <> 4)) then
       dwPdr:=dwPdr+ R[Zu[l]]*fi(1,ys_real[k,l],-ys_imag[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
     if ((l <> k) and (mhu[l] <> 0)) then
       dwQdr:=dwQdr+ R[Zf[l]]*ytrue[l]*fi(3,ys_real[k,l],-ys_imag[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
    end;
   end;
  end
 else
  begin
   if (mhu[k] <> 4) then
     dwPdr:=dwPdr+ R[Zu[k]]*fi(1,ys_real[k,i],-ys_imag[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
   dwQdr:=dwQdr+ R[Zf[k]]*ytrue[k]*fi(2,ys_real[k,i],-ys_imag[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
   dwQdr:=dwQdr+ R[Zf[i]]*ytrue[k]*fi(3,ys_real[k,i],-ys_imag[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
  end;
 dwPdr :=dwPdr*Kn;
 dwQdr :=dwQdr*Kn;
 aa[kJac+(n div 2),iJac] := dwPdr+dwQdr;
 iJac:=iJac+1;
 end;
 kJac:=kJac+1;
 end;
iJac:=1;
kJac:=1;


//dVq/dU
for k:=1 to ntrue do
 if ((mhu[k] <> 0) and (mhu[k] <> 4)) then
 begin
 iJac:=1;
 for i:=1 to ntrue do
 if ((mhu[i] <> 0) and (mhu[i] <> 4)) then
 begin
 dwPdr := 0;
 dwQdr := 0;
 if (i=k) then
  begin
   dwPdr:=dwPdr+ 2.*ys_imag[k,k]*R[Zu[k]];
   for j:=1 to ntrue do
    if (j <> k) then
     dwQdr:=dwQdr+R[Zf[k]]*ytrue[j]*fi(2,ys_imag[k,j],ys_real[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
   for l:=1 to ntrue do
    if ((l <> k) and (mhu[l] <> 0)) then
     begin
      dwQdr:=dwQdr+ R[Zf[l]]*ytrue[l]*fi(3,ys_imag[k,l],ys_real[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
      if (mhu[l] <> 4) then
       dwPdr:=dwPdr+ R[Zu[l]]*fi(1,ys_imag[k,l],ys_real[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
     end;
  end
 else
  begin
   dwPdr:=dwPdr+ R[Zu[k]]*fi(1,ys_imag[k,i],ys_real[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
   dwQdr:=dwQdr+ R[Zf[k]]*ytrue[k]*fi(2,ys_imag[k,i],ys_real[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
   dwQdr:=dwQdr+ R[Zf[i]]*ytrue[k]*fi(3,ys_imag[k,i],ys_real[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
  end;
 dwPdr :=dwPdr* Kn;
 dwQdr :=dwQdr* Kn;
 aa[kJac+(n div 2)+ntrue-nb,iJac] := dwPdr+dwQdr;
 iJac:=iJac+1;
 end;
 kJac:=kJac+1;
end;
//***************************
kJac:=1;
iJac:=1;

//���� ��� dP/dD
for k :=1 to ntrue do
 if (mhu[k] <> 0) then
 begin
 iJac:=1;
 for i :=1 to ntrue do
 if (mhu[i] <> 0) then
 begin
  dwPdi :=0;
  if (i=k) then
  begin
  for j :=1 to ntrue do
    begin
        if (j <> k) then
         begin
          dwPdi:=dwPdi+ ytrue[k]*ytrue[j]*(ys_real[k,j]*(-sin(ytrue[k+ntrue])*cos(ytrue[j+ntrue])+cos(ytrue[k+ntrue])*sin(ytrue[j+ntrue]))-ys_imag[k,j]*(-sin(ytrue[k+ntrue])*sin(ytrue[j+ntrue])-cos(ytrue[k+ntrue])*cos(ytrue[j+ntrue])));
         end;
    end;
  end
  else
  begin
   dwPdi:=dwPdi+ ytrue[k]*ytrue[i]*(ys_real[k,i]*(-cos(ytrue[k+ntrue])*sin(ytrue[i+ntrue])+sin(ytrue[k+ntrue])*cos(ytrue[i+ntrue]))-ys_imag[k,i]*(cos(ytrue[k+ntrue])*cos(ytrue[i+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[i+ntrue])));
  end;
 dwPdi:=dwPdi* Kn;
 aa[kJac,iJac+ntrue-nb-ng] := dwPdi;
 iJac:=iJac+1;
 end;
 kJac:=kJac+1;
 end;
iJac:=1;
kJac:=1;


//���� ���  dQ/dD
for k :=1 to ntrue do
 if ((mhu[k] <> 0) and (mhu[k] <> 4)) then
 begin
 iJac:=1;
 for i :=1 to ntrue do
 if (mhu[i] <> 0) then
 begin
  dwQdi :=0;
  if (i=k) then
  begin
  for j :=1 to ntrue do
    begin
        if (j <> k) then
         begin
          dwQdi:=dwQdi+ ytrue[k]*ytrue[j]*(ys_imag[k,j]*(-sin(ytrue[k+ntrue])*cos(ytrue[j+ntrue])+cos(ytrue[k+ntrue])*sin(ytrue[j+ntrue]))+ys_real[k,j]*(-sin(ytrue[k+ntrue])*sin(ytrue[j+ntrue])-cos(ytrue[k+ntrue])*cos(ytrue[j+ntrue])));
         end;
    end;
  end
  else
  begin
   dwQdi:=dwQdi+ ytrue[k]*ytrue[i]*(ys_imag[k,i]*(-cos(ytrue[k+ntrue])*sin(ytrue[i+ntrue])+sin(ytrue[k+ntrue])*cos(ytrue[i+ntrue]))+ys_real[k,i]*(cos(ytrue[k+ntrue])*cos(ytrue[i+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[i+ntrue])));
  end;
 dwQdi:=dwQdi* Kn;
 aa[kJac+ntrue-nb,iJac+ntrue-nb-ng] := dwQdi;
 iJac:=iJac+1;
 end;
 kJac:=kJac+1;
 end;

kJac:=1;
iJac:=1;
//***********************
//���� ��� dVp/df
for k:=1 to ntrue do
 if (mhu[k] <> 0) then
 begin
 iJac:=1;
 for i:=1 to ntrue do
 if (mhu[i] <> 0) then
 begin
 dwPdr := 0;
 dwQdr := 0;
 if (i=k) then
  begin
   if (mhu[k] <> 4) then
   begin
    for j:=1 to ntrue do
     begin
      if (j <> k) then
       begin
        dwPdr:=dwPdr+ R[Zu[k]]*ytrue[j]*fi(2,ys_real[k,j],-ys_imag[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
        dwQdr:=dwQdr+ R[Zf[k]]*ytrue[k]*ytrue[j]*fi(4,ys_real[k,j],-ys_imag[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
       end;
     end;
    for l:=1 to ntrue do
     begin
      if ((l <> k) and (mhu[l] <> 0)) then
       begin
        dwQdr:=dwQdr+ R[Zf[l]]*ytrue[k]*ytrue[l]*fi(1,ys_real[k,l],-ys_imag[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
        if (mhu[l] <> 4) then
           dwPdr:=dwPdr+ R[Zu[l]]*ytrue[k]*fi(2,ys_real[k,l],-ys_imag[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
       end;
     end;
   end
   else //��� ������������ �����
   begin
    for j:=1 to ntrue do
     begin
      if (j <> k) then
        dwQdr:=dwQdr+ R[Zf[k]]*ytrue[k]*ytrue[j]*fi(4,ys_real[k,j],-ys_imag[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
     end;
    for l:=1 to ntrue do
     begin
      if ((mhu[l] <> 0) and (mhu[l] <> 4)) then
        dwPdr:=dwPdr+ R[Zu[l]]*ytrue[k]*fi(2,ys_real[k,l],-ys_imag[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
      if ((l <> k) and (mhu[l] <> 0)) then
        dwQdr:=dwQdr+ R[Zf[l]]*ytrue[k]*ytrue[l]*fi(1,ys_real[k,l],-ys_imag[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
     end;
   end;
  end
 else
  begin
   if (mhu[k] <> 4) then
     dwPdr:=dwPdr+ R[Zu[k]]*ytrue[i]*fi(3,ys_real[k,i],-ys_imag[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
   dwQdr:=dwQdr+ R[Zf[k]]*ytrue[k]*ytrue[i]*fi(1,ys_real[k,i],-ys_imag[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
   dwPdr:=dwPdr+ R[Zu[i]]*ytrue[k]*fi(3,ys_real[k,i],-ys_imag[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
   dwQdr:=dwQdr+ R[Zf[i]]*ytrue[k]*ytrue[i]*fi(4,ys_real[k,i],-ys_imag[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
  end;
 dwPdr :=dwPdr* Kn;
 dwQdr :=dwQdr* Kn;
 aa[kJac+(n div 2),iJac+ntrue-nb-ng] := dwPdr+dwQdr;
 iJac:=iJac+1;
 end;
 kJac:=kJac+1;
end;

iJac:=1;
kJac:=1;


//dVq/df
for k:=1 to ntrue do
 if ((mhu[k] <> 0) and (mhu[k] <> 4)) then
 begin
 iJac:=1;
 for i:=1 to ntrue do
 if (mhu[i] <> 0) then
 begin
 dwPdr := 0;
 dwQdr := 0;
 if (i=k) then
  begin
   for j:=1 to ntrue do
    begin
     if (j <> k) then
      begin
       dwPdr:=dwPdr+R[Zu[k]]*ytrue[j]*fi(2,ys_imag[k,j],ys_real[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
       dwQdr:=dwQdr+ R[Zf[k]]*ytrue[k]*ytrue[j]*fi(4,ys_imag[k,j],ys_real[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
      end;
    end;
   for l:=1 to ntrue do
    begin
     if((l <> k) and (mhu[l] <> 0)) then
      begin
       dwQdr:=dwQdr+ R[Zf[l]]*ytrue[k]*ytrue[l]*fi(1,ys_imag[k,l],ys_real[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
       if (mhu[l] <> 4) then
          dwPdr:=dwPdr+ R[Zu[l]]*ytrue[k]*fi(2,ys_imag[k,l],ys_real[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
      end;
    end;
  end
 else
  begin
   dwPdr:=dwPdr+ R[Zu[k]]*ytrue[i]*fi(3,ys_imag[k,i],ys_real[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
   dwQdr:=dwQdr+ R[Zf[k]]*ytrue[k]*ytrue[i]*fi(1,ys_imag[k,i],ys_real[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
   dwPdr:=dwPdr+ R[Zu[i]]*ytrue[k]*fi(3,ys_imag[k,i],ys_real[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
   dwQdr:=dwQdr+ R[Zf[i]]*ytrue[k]*ytrue[i]*fi(4,ys_imag[k,i],ys_real[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
  end;
 dwPdr := dwPdr*Kn;
 dwQdr := dwQdr*Kn;
 aa[kJac+(n div 2)+ntrue-nb,iJac+ntrue-nb-ng] := dwPdr+dwQdr;
 iJac:=iJac+1;
 end;
 kJac:=kJac+1;
end;
//***************************
kJac:=1;
iJac:=1;

//dVp/dRu
for k:=1 to ntrue do
 if (mhu[k] <> 0) then
 begin
 iJac:=1;
 for i:=1 to ntrue do
 if ((mhu[i] <> 0) and (mhu[i] <> 4)) then
 begin
 dwPdr := 0;
 if (k=i) then
  begin
   dwPdr:=dwPdr+ 2.*ys_real[k,k]*ytrue[k];
   for j:=1 to ntrue do
    begin
     if (j <> k) then
       dwPdr:=dwPdr+ytrue[j]*fi(1,ys_real[k,j],-ys_imag[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
    end;
  end
 else
 begin
  dwPdr:=dwPdr+ ytrue[k]*fi(1,ys_real[k,i],-ys_imag[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
 end;
 dwPdr :=dwPdr* Kn;
 aa[kJac+(n div 2),iJac+(n div 2)] := dwPdr;
 iJac:=iJac+1;
 end;
 kJac:=kJac+1;
 end;

kJac:=1;
iJac:=1;

//dVQ/dRu
for k:=1 to ntrue do
 if ((mhu[k] <> 0) and (mhu[k] <> 4)) then
 begin
 iJac:=1;
 for i:=1 to ntrue do
 if ((mhu[i] <> 0) and (mhu[i] <> 4)) then
 begin
 dwQdr := 0;
 if (k=i) then
  begin
   dwQdr:=dwQdr+ 2.*ys_imag[k,k]*ytrue[k];
   for j:=1 to ntrue do
    begin
     if (j <> k) then
       dwQdr:=dwQdr+ ytrue[j]*fi(1,ys_imag[k,j],ys_real[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
    end;
  end
 else
 begin
  dwQdr:=dwQdr+ ytrue[k]*fi(1,ys_imag[k,i],ys_real[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
 end;
 dwQdr :=dwQdr* Kn;
 aa[kJac+(n div 2)+ntrue-nb,iJac+(n div 2)] := dwQdr;
 iJac:=iJac+1;
 end;
 kJac:=kJac+1;
 end;


iJac:=1;
kJac:=1;
//dVp/dRf
for k:=1 to ntrue do
 if (mhu[k] <> 0) then
 begin
 iJac:=1;
 for i:=1 to ntrue do
 if (mhu[i] <> 0) then
 begin
 dwPdr := 0;
 if (k=i) then
  begin
   for j:=1 to ntrue do
    begin
     if (j <> k) then
       dwPdr:=dwPdr+ ytrue[k]*ytrue[j]*fi(2,ys_real[k,j],-ys_imag[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
    end;
  end
 else
 begin
  dwPdr:=dwPdr+ ytrue[k]*ytrue[i]*fi(3,ys_real[k,i],-ys_imag[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
 end;
 dwPdr :=dwPdr* Kn;
 aa[kJac+(n div 2),iJac+(n div 2)+ntrue-nb-ng] := dwPdr;
 iJac:=iJac+1;
 end;
 kJac:=kJac+1;
 end;

kJac:=1;
iJac:=1;

//dVq/dRf
for k:=1 to ntrue do
 if ((mhu[k] <> 0) and (mhu[k] <> 4)) then
 begin
 iJac:=1;
 for i:=1 to ntrue do
 if (mhu[i] <> 0) then
 begin
 dwQdr := 0;
 if (k=i) then
  begin
   for j:=1 to ntrue do
    begin
     if (j <> k) then
       dwQdr:=dwQdr+ ytrue[k]*ytrue[j]*fi(2,ys_imag[k,j],ys_real[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
    end;
  end
 else
 begin
  dwQdr:=dwQdr+ ytrue[k]*ytrue[i]*fi(3,ys_imag[k,i],ys_real[k,i],ytrue[ntrue+k],ytrue[ntrue+i]);
 end;
 dwQdr :=dwQdr*Kn;
 aa[kJac+(n div 2)+ntrue-nb,iJac+(n div 2)+ntrue-nb-ng] := dwQdr;
 iJac:=iJac+1;
 end;
 kJac:=kJac+1;
 end;

//dfp/dKn
kJac := 1;
for k := 1 to ntrue do
 if (mhu[k] <> 0) then
 begin
//  Nebal[kJac]=0;
  tmp1 := 0;
  tmp2 := 0;
  tmp1 := -sng_real[k];
  tmp2 :=tmp2+ ys_real[k,k]*ytrue[k]*ytrue[k];
  for j := 1 to ntrue do
    if (j <>k) then
    begin
     tmp2 :=tmp2+ ytrue[k]*ytrue[j]*(ys_real[k,j]*(cos(ytrue[ntrue+k])*cos(ytrue[ntrue+j])+sin(ytrue[ntrue+k])*sin(ytrue[ntrue+j]))-ys_imag[k,j]*(cos(ytrue[ntrue+k])*sin(ytrue[ntrue+j])-sin(ytrue[ntrue+k])*cos(ytrue[ntrue+j])));
    end;
  aa[kJac,n] := tmp2;
//  eql[kJac] := tmp1+tmp2*Kn;
 // Nebal[kJac]=tmp1+tmp2;
  kJac:=kJac+1;
 end;

//���� ��� Q
kJac := 1;
for k := 1 to ntrue do
if ((mhu[k] <> 0) and (mhu[k] <> 4)) then
 begin
//  Nebal[kJac+n-nb]=0;
  tmp2 := 0;
  tmp1 := 0;
  tmp1 := -sng_imag[k];
  tmp2 :=tmp2+ ys_imag[k,k]*ytrue[k]*ytrue[k];
  for j := 1 to ntrue do
    if (j <> k) then
    begin
    tmp2 :=tmp2+ ytrue[k]*ytrue[j]*(ys_imag[k,j]*(cos(ytrue[ntrue+k])*cos(ytrue[ntrue+j])+sin(ytrue[ntrue+k])*sin(ytrue[ntrue+j]))+ys_real[k,j]*(cos(ytrue[ntrue+k])*sin(ytrue[ntrue+j])-sin(ytrue[ntrue+k])*cos(ytrue[ntrue+j])));
    end;
  aa[kJac+ntrue-nb,n] := tmp2;
//  eql[kJac+n-nb] := tmp1+tmp2*Kn;
//  Nebal[kJac+n-nb]=tmp1+tmp2;
  kJac:=kJac+1;
 end;

//���� ��� Vp
kJac := 1;
for k := 1 to ntrue do
 if (mhu[k] <> 0) then
 begin
  tmp1 := 0;
  tmp2 := 0;
  for l := 1 to ntrue do
   if((mhu[l] <> 0) and (mhu[l] <> 4)) then
   begin
    if (k=l) then
     begin
      tmp1 :=tmp1+ 2.*ys_real[k,k]*ytrue[k]*R[Zu[k]];
      for j := 1 to ntrue do
       if (j <> k) then
         tmp1 :=tmp1+ R[Zu[l]]*ytrue[j]*fi(1,ys_real[k,j],-ys_imag[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
     end
    else
      tmp1 :=tmp1+ R[Zu[l]]*ytrue[k]*fi(1,ys_real[k,l],-ys_imag[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
   end;
  for l := 1 to ntrue do
   if (mhu[l] <> 0) then
   begin
    if (k=l) then
     begin
      for j := 1 to ntrue do
       if (j <> k) then
         tmp2 :=tmp2+ R[Zf[l]]*ytrue[k]*ytrue[j]*fi(2,ys_real[k,j],-ys_imag[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
     end
    else
      tmp2 :=tmp2+ R[Zf[l]]*ytrue[k]*ytrue[l]*fi(3,ys_real[k,l],-ys_imag[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
   end;
 aa[(n div 2)+kJac,n] := tmp1+tmp2;
// eql[(n div 2)+kJac] := (tmp1*Kn+tmp2*Kn);
 kJac:=kJac+1;
 end;

//Vq
kJac := 1;
for k := 1 to ntrue do
 if ((mhu[k] <> 0)and (mhu[k]<>4)) then
 begin
  tmp1 := 0;
  tmp2 := 0;
  for l := 1 to ntrue do
   if((mhu[l] <> 0) and (mhu[l] <> 4)) then
   begin
    if (k=l) then
     begin
      tmp1 :=tmp1+ 2.*ys_imag[k,k]*ytrue[k]*R[Zu[k]];
      for j := 1 to ntrue do
       if (j <> k) then
         tmp1 :=tmp1+ R[Zu[l]]*ytrue[j]*fi(1,ys_imag[k,j],ys_real[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
     end
    else
      tmp1 :=tmp1+ R[Zu[l]]*ytrue[k]*fi(1,ys_imag[k,l],ys_real[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
   end;
  for l := 1 to ntrue do
   if (mhu[l] <> 0) then
   begin
    if (k=l) then
     begin
      for j := 1 to ntrue do
       if (j <> k) then
         tmp2 :=tmp2+ R[Zf[l]]*ytrue[k]*ytrue[j]*fi(2,ys_imag[k,j],ys_real[k,j],ytrue[ntrue+k],ytrue[ntrue+j]);
     end
    else
      tmp2 :=tmp2+ R[Zf[l]]*ytrue[k]*ytrue[l]*fi(3,ys_imag[k,l],ys_real[k,l],ytrue[ntrue+k],ytrue[ntrue+l]);
   end;
 aa[(n div 2)+kJac+ntrue-nb,n] := tmp1+tmp2;
// eql[(n div 2)+kJac+ntrue-nb] := (tmp1*Kn+tmp2*Kn);
 kJac:=kJac+1;
 end;

{


kJac:=1;
for k:=1 to ntrue do
begin
 if(mhu[k]<>0)then
 begin
 iJac:=1;
 for i:=1 to ntrue do
 if((mhu[i]<>0)and(mhu[i]<>4))then
  begin
  dwPdr:=0;
  if(i=k)then
  begin
  dwPdr:=dwPdr+2.*ys_real[k,k]*ytrue[k];
  for j:=1 to ntrue do
    begin
     if(j<>k)then
       dwPdr:=dwPdr+ytrue[j]*(ys_real[k,j]*(cos(ytrue[k+ntrue])*cos(ytrue[j+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[j+ntrue]))-ys_imag[k,j]*(cos(ytrue[k+ntrue])*sin(ytrue[j+ntrue])-sin(ytrue[k+ntrue])*cos(ytrue[j+ntrue])));
    end;
  end
  else
   begin
    dwPdr:=dwPdr+ytrue[k]*(ys_real[k,i]*(cos(ytrue[k+ntrue])*cos(ytrue[i+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[i+ntrue]))-ys_imag[k,i]*(cos(ytrue[k+ntrue])*sin(ytrue[i+ntrue])-sin(ytrue[k+ntrue])*cos(ytrue[i+ntrue])));
   end;
  aa[kJac,iJac]:=dwPdr;
  iJac:=iJac+1;
  end;
 kJac:=kJac+1;
 end;
end;
iJac:=1;
kJac:=1;

//���� ��� dQ/dU
for k:=1 to ntrue do
begin
 if((mhu[k]<>0)and(mhu[k]<>4))then
 begin
 iJac:=1;
 for i:=1 to ntrue do
 if((mhu[i]<>0)and(mhu[i]<>4))then
  begin
  dwQdr:=0;
  if(i=k)then
  begin
  dwQdr:=dwQdr+2.*ys_imag[k,k]*ytrue[k];
  for j:=1 to ntrue do
    begin
     if(j<>k)then
       dwQdr:=dwQdr+ytrue[j]*(ys_imag[k,j]*(cos(ytrue[k+ntrue])*cos(ytrue[j+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[j+ntrue]))+ys_real[k,j]*(cos(ytrue[k+ntrue])*sin(ytrue[j+ntrue])-sin(ytrue[k+ntrue])*cos(ytrue[j+ntrue])));
    end;
  end
  else
   begin
    dwQdr:=dwQdr+ytrue[k]*(ys_imag[k,i]*(cos(ytrue[k+ntrue])*cos(ytrue[i+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[i+ntrue]))+ys_real[k,i]*(cos(ytrue[k+ntrue])*sin(ytrue[i+ntrue])-sin(ytrue[k+ntrue])*cos(ytrue[i+ntrue])));
   end;
  aa[kJac+ntrue-nb,iJac]:=dwQdr;
  iJac:=iJac+1;
  end;
 kJac:=kJac+1;
 end;
end;


//���� ��� dP/dD
kJac:=1;
for k:=1 to ntrue do
begin
 if(mhu[k]<>0)then
 begin
 iJac:=1;
 for i:=1 to ntrue do
 if(mhu[i]<>0)then
  begin
  dwPdi:=0;
  if(i=k)then
  begin
  for j:=1 to ntrue do
    begin
     if(j<>k)then
       dwPdi:=dwPdi+ytrue[k]*ytrue[j]*(ys_real[k,j]*(-sin(ytrue[k+ntrue])*cos(ytrue[j+ntrue])+cos(ytrue[k+ntrue])*sin(ytrue[j+ntrue]))-ys_imag[k,j]*(-sin(ytrue[k+ntrue])*sin(ytrue[j+ntrue])-cos(ytrue[k+ntrue])*cos(ytrue[j+ntrue])));
    end;
  end
  else
   begin
    dwPdi:=dwPdi+ytrue[k]*ytrue[i]*(ys_real[k,i]*(-cos(ytrue[k+ntrue])*sin(ytrue[i+ntrue])+sin(ytrue[k+ntrue])*cos(ytrue[i+ntrue]))-ys_imag[k,i]*(cos(ytrue[k+ntrue])*cos(ytrue[i+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[i+ntrue])));
   end;
  aa[kJac,iJac+ntrue-nb-ng]:=dwPdi;
  iJac:=iJac+1;
  end;
 kJac:=kJac+1;
 end;
end;


//���� ���  dQ/dD
kJac:=1;
for k:=1 to ntrue do
begin
 if((mhu[k]<>0)and(mhu[k]<>4))then
 begin
 iJac:=1;
 for i:=1 to ntrue do
 if(mhu[i]<>0)then
  begin
  dwQdi:=0;
  if(i=k)then
  begin
  for j:=1 to ntrue do
    begin
     if(j<>k)then
       dwQdi:=dwQdi+ytrue[k]*ytrue[j]*(ys_imag[k,j]*(-sin(ytrue[k+ntrue])*cos(ytrue[j+ntrue])+cos(ytrue[k+ntrue])*sin(ytrue[j+ntrue]))+ys_real[k,j]*(-sin(ytrue[k+ntrue])*sin(ytrue[j+ntrue])-cos(ytrue[k+ntrue])*cos(ytrue[j+ntrue])));
    end;
  end
  else
   begin
    dwQdi:=dwQdi+ytrue[k]*ytrue[i]*(ys_imag[k,i]*(-cos(ytrue[k+ntrue])*sin(ytrue[i+ntrue])+sin(ytrue[k+ntrue])*cos(ytrue[i+ntrue]))+ys_real[k,i]*(cos(ytrue[k+ntrue])*cos(ytrue[i+ntrue])+sin(ytrue[k+ntrue])*sin(ytrue[i+ntrue])));
   end;
  aa[kJac+ntrue-nb,iJac+ntrue-nb-ng]:=dwQdi;
  iJac:=iJac+1;
  end;
 kJac:=kJac+1;
 end;
end;
}
end;

procedure nldet(n:integer;var a:nla;var d:real);
  label E;
  var max,aa:real;
  var i,j,jj,l,k,ko:integer;
        begin

  if n=1 then begin d:=a[n,n];goto E end;
  ko:=0;d:=1.;
  for i:=1 to n-1 do
        begin
  k:=i;l:=i;max:=a[i,i];
  for j:=i to n do
        begin
  for jj:=i to n do
        begin
  if(abs(max)-abs(a[j,jj]))<0. then
  begin max:=a[j,jj];k:=j;l:=jj end
         end
        end;
  if max=0. then begin d:=0.;goto E end;
  if k<>i then
        begin
   ko:=ko+1;
   for j:=1 to n do
   begin aa:=a[k,j];a[k,j]:=a[i,j];a[i,j]:=aa end;
        end;
   if l<>i then
        begin
   ko:=ko+1;
   for j:=1 to n do
   begin aa:=a[j,l];a[j,l]:=a[j,i];a[j,i]:=aa end;
        end;
   for j:=i+1 to n do
        begin
   aa:=a[j,i]/a[i,i];
   for jj:=i to n do a[j,jj]:=a[j,jj]-a[i,jj]*aa;
         end;
        end;
  for i:=1 to n do d:=d*a[i,i];
  if(ko mod 2)>0 then d:=-d;
  E:    end;

procedure nleqc(ndim:integer;var x,direk:real;var y,yy,y0,{eql,}dery,prmt:nly;var aa,aap:nla);
    var k1,k2,k0:integer;
    var s,d:real;
    begin
    {f1(eql,y0);}
    f2(aa,y);
    s:=0.;
    for k1:=1 to ndim-1 do
    begin for k2:=1 to ndim-1 do aap[k1,k2]:=aa[k1,k2] end;
    for k0:=1 to ndim-1 do
         begin
    for k2:=1 to ndim-1 do aa[k2,k0]:=eql[k2];
    nldet(ndim-1,aa,d); dery[k0]:=direk*d;
    for k1:=1 to ndim-1 do
    begin for k2:=1 to ndim-1 do aa[k1,k2]:=aap[k1,k2] end;
         end;
    nldet(ndim-1,aa,d); dery[ndim]:=direk*d;
    for k1:=1 to ndim do s:=s+abs(dery[k1]);
    if s=0. then begin s:=1.; for k1:=1 to ndim-1 do y[k1]:=y[k1]+0.00000001;end;
    for k1:=1 to ndim do dery[k1]:=dery[k1]/s;
end;

procedure nleqd(ndim,ihlf:integer;var x,direk:real;var y,yy,{eql,}prmt:nly;var aa,aap:nla);
    label w;
    var    i:integer;
           s:real;

begin
   if ihlf>=10 then
    MessageDlg('Reduce the requirement to admissible error', mtError, [mbOk], 0);
    s:=y[ndim]*yy[ndim];
  if s>0. then goto w;
    {f1({eql,y);}
    jk:=jk+1;
    if (jk>korn) and gg then
      begin
       gg:=false;
       messagedlg('More than 200 approximations were detected',mtInformation, [mbOk], 0)
      end;
    if jk<=korn then
    begin
    for i:=1 to n do
     begin
     aprl[i,jk]:= y[i];
     aprr[i,jk]:=yy[i];
     Kn:=Kn-1+1;
     end;
    end;
w:   yy:=y;
end;

procedure nlmain(ndim:integer;var ihlf:integer;var y0,dery,eql,y,yy,prmt:nly;
                  var aa,aap:nla;var aux:nlaux;var direk:real);
    label 4,9,11,13,15,21,23,29,100,200,201,206,210,222,111;
    var delt,x,z,h:real;
    var istep,isw,n,i:integer;
 begin
   n:=1;ihlf:=0;x:=prmt[1];h:=prmt[3];prmt[5]:=0.;y[ndim]:=1.;yy[ndim]:=y[ndim];
   for i:=1 to ndim do
        begin aux[16,i]:=0.; aux[15,i]:=dery[i];aux[1,i]:=y[i] end;
   if (h*(prmt[2]-x))<0. then ihlf:=13;
   if (h*(prmt[2]-x))=0. then ihlf:=12;
4: nleqc(ndim,x,direk,y,yy,y0,{eql,}dery,prmt,aa,aap);
   nleqd(ndim,ihlf,x,direk,y,yy,{eql,}prmt,aa,aap);
   if(prmt[5]<0.0) or (prmt[5]>0.0) then goto 111;
   if(prmt[5]=0.0) and (ihlf>0)    then goto 111;
   for i:=1 to ndim do aux[8,i]:=dery[i];
   isw:=1;
   goto 100;
9: x:=x+h;
   for i:=1 to ndim do aux[2,i]:=y[i];
11: ihlf:=ihlf+1; x:=x-h;
    for i:=1 to ndim do aux[4,i]:=aux[2,i];
    h:=0.5*h; n:=1;
    isw:=2;
    goto 100;
13: x:=x+h;
    nleqc(ndim,x,direk,y,yy,y0,{eql,}dery,prmt,aa,aap);
    n:=2;
    for i:=1 to ndim do begin aux[2,i]:=y[i]; aux[9,i]:=dery[i] end;
    isw:=3;
    goto 100;
15: delt:=0.;
    for i:=1 to ndim do delt:=delt+aux[15,i]*abs(y[i]-aux[4,i]);
    delt:=0.06666667*delt;
    if(delt-prmt[4])>0. then
        begin
    if(ihlf-10)<0 then goto 11;
    ihlf:=11; x:=x+h;
    goto 4
        end;x:=x+h;
    nleqc(ndim,x,direk,y,yy,y0,{eql,}dery,prmt,aa,aap);
    for i:=1 to ndim do begin aux[3,i]:=y[i]; aux[10,i]:=dery[i] end;
    n:=3;
    isw:=4;
    goto 100;
21: n:=1; x:=x+h;
    nleqc(ndim,x,direk,y,yy,y0,{eql,}dery,prmt,aa,aap);
    x:=prmt[1];
    for i:=1 to ndim do
         begin
    aux[11,i]:=dery[i];
    y[i]:=aux[1,i]+h*(0.375*aux[8,i]+0.7916667*aux[9,i]-0.2083333*aux[10,i]+0.04166667*dery[i])
         end;
23: x:=x+h; n:=n+1;
    nleqc(ndim,x,direk,y,yy,y0,{eql,}dery,prmt,aa,aap);
    nleqd(ndim,ihlf,x,direk,y,yy,{eql,}prmt,aa,aap);
    if(prmt[5]<0.0) or (prmt[5]<0.0) then goto 111;
    if(n-4)>=0 then goto 200;
    for i:=1 to ndim do begin aux[n,i]:=y[i]; aux[n+7,i]:=dery[i] end;
    if(n-3)>0 then goto 200;
    if(n-3)=0 then goto 29;
    for i:=1 to ndim do
         begin
    delt:=aux[9,i]+aux[9,i]; delt:=delt+delt;
    y[i]:=aux[1,i]+0.3333333*h*(aux[10,i])
         end;
    goto 23;
29: for i:=1 to ndim do begin delt:=aux[9,i]+aux[10,i]; delt:=delt+delt+delt;
                              y[i]:=aux[1,i]+0.375*h*(aux[8,i]+delt+aux[11,i])
                        end;
    goto 23;
100:for i:=1 to ndim do begin z:=h*aux[n+7,i];aux[5,i]:=z;
                              y[i]:=aux[n,i]+0.4*z;
                        end;
    z:=x+0.4*h;
    nleqc(ndim,z,direk,y,yy,y0,{eql,}dery,prmt,aa,aap);
    for i:=1 to ndim do begin z:=h*dery[i]; aux[6,i]:=z;
                              y[i]:=aux[n,i]+0.2969776*aux[5,i]+0.1587596*z
                        end;
    z:=x+0.4557372*h;
    nleqc(ndim,z,direk,y,yy,y0,{eql,}dery,prmt,aa,aap);
    for i:=1 to ndim do
    begin z:=h*dery[i];aux[7,i]:=z;
          y[i]:=aux[n,i]+0.2181004*aux[5,i]-3.050965*aux[6,i]*3.832865*z
    end;
    z:=x+h;
    nleqc(ndim,z,direk,y,yy,y0,{eql,}dery,prmt,aa,aap);
    for i:=1 to ndim do
    y[i]:=aux[n,i]+0.1747603*aux[5,i]-0.5514807*aux[6,i]+1.205536*aux[7,i]+0.17118448*h*dery[i];
    if isw=1 then goto 9;
    if isw=2 then goto 13;
    if isw=3 then goto 15;
    if isw=4 then goto 21;
200:istep:=3;
201:if (n-8)=0 then
       begin
    for n:=2 to 7 do
       begin
    for i:=1 to ndim do begin aux[n-1,i]:=aux[n,i];aux[n+6,i]:=aux[n+7,i] end;
       end;
    n:=7
       end;
    n:=n+1;
    for i:=1 to ndim do begin aux[n-1,i]:=y[i]; aux[n+6,i]:=dery[i] end;
    x:=x+h;
206:istep:=istep+1;
    for i:=1 to ndim do
       begin
    delt:=aux[n-4,i]+1.333333*h*(aux[n+6,i]+aux[n+6,i]-aux[n+5,i]+aux[n+4,i]+aux[n+4,i]);
    y[i]:=delt-0.9256198*aux[16,i];
    aux[16,i]:=delt
       end;
    nleqc(ndim,x,direk,y,yy,y0,{eql,}dery,prmt,aa,aap);
            for i:=1 to ndim do
       begin
    delt:=0.125*(9.*aux[n-1,i]-aux[n-3,i]+3.*h*(dery[i]+aux[n+6,i]+aux[n+6,i]-aux[n+5,i]));
    aux[16,i]:=aux[16,i]-delt;
    y[i]:=delt+0.07438017*aux[16,i]
       end;
    delt:=0.;
    for i:=1 to ndim do delt:=delt+aux[15,i]*abs(aux[16,i]);
    if (delt-prmt[4])>=0. then goto 222;
210:
    nleqc(ndim,x,direk,y,yy,y0,{eql,}dery,prmt,aa,aap);
    nleqd(ndim,ihlf,x,direk,y,yy,{eql,}prmt,aa,aap);
    if prmt[5]<>0. then goto 111;
    if(ihlf-11)>=0 then goto 111;
    if(h*(x-prmt[2]))>=0. then goto 111;
    if(abs(x-prmt[2])-0.1*abs(h))<0. then goto 111;
    if(delt-0.02*prmt[4])>0. then goto 201;
    if ihlf<=0 then goto 201;
    if (n-7)<0 then goto 201;
    if (istep-4)<0 then goto 201;
    if(istep mod 2 )<>0 then goto 201;
    h:=h+h;
    ihlf:=ihlf-1;
    istep:=0;
    for i:=1 to ndim do
         begin
    aux[n-1,i]:=aux[n-2,i];aux[n-2,i]:=aux[n-4,i];
    aux[n-3,i]:=aux[n-6,i];aux[n+6,i]:=aux[n+5,i];
    aux[n+5,i]:=aux[n+3,i];aux[n+4,i]:=aux[n+1,i];
    delt:=aux[n+6,i]+aux[n+5,i];
    delt:=delt+delt+delt;
    aux[16,i]:=8.962963*(y[i]-aux[n-3,i])-3.361111*h*(dery[i]+delt+aux[n+4,i])
         end;
    goto 201;
222:ihlf:=ihlf+1;
    if(ihlf-10)>0 then goto 210;
    h:=0.5*h;
    istep:=0;
    for i:=1 to ndim do
         begin
    y[i]:=0.00390625*(80.*aux[n-1,i]+135.*aux[n-2,i]+40.*aux[n-3,i]+
          aux[n-4,i])-0.1171875*(aux[n+6,i]-6.*aux[n+5,i]-aux[n+4,i])*h;
    aux[n-4,i]:=0.00390625*(12.*aux[n-1,i]+135.*aux[n-2,i]+
          108.*aux[n-3,i]+aux[n-4,i])-0.0234375*(aux[n+6,i]+18.*aux[n+5,i]-9.*aux[n+4,i])*h;
    aux[n-3,i]:=aux[n-2,i];
    aux[n+4,i]:=aux[n+5,i]
          end;
    x:=x-h;
    delt:=x-(h+h);
    nleqc(ndim,delt,direk,y,yy,y0,{eql,}dery,prmt,aa,aap);
    for i:=1 to ndim do
          begin
    aux[n-2,i]:=y[i];
    aux[n+5,i]:=dery[i];
    y[i]:=aux[n-4,i]
          end;
    delt:=delt-(h+h);
    nleqc(ndim,delt,direk,y,yy,y0,{eql,}dery,prmt,aa,aap);
    for i:=1 to ndim do
          begin
    delt:=aux[n+5,i]+aux[n+4,i];
    delt:=delt+delt+delt;
    aux[16,i]:=8.962963*(aux[n-1,i]-y[i])-3.361111*h*(aux[n+6,i]+delt+dery[i]);
    aux[n+3,i]:=dery[i]
          end;
    goto 206;
111:
          end;

procedure dragmain;
var i,ii:integer;
begin
  f1({eql,}y0);
  prmt[1]:=0.;
  for ii:=1 to 2 do
     begin direk:=-1.; if ii=2 then direk:=1.;
  for i:=1 to n+1 do
    begin dery[i]:=1./(n+1);y[i]:=y0[i];yy[i]:=y[i] end;
    nlmain(n+1,ihlf,y0,dery,eql,y,yy,prmt,a,aap,aux,direk);
 end;
end;


function fi(i:integer;y1:real;y2:real;f1:real;f2:real):real;
begin
case (i) of
1:
  begin
   (* C2PAS: Exit *) Result := y1*(cos(f1)*cos(f2)+sin(f1)*sin(f2))+y2*(cos(f1)*sin(f2)-sin(f1)*cos(f2));
  end;
2:
  begin
   (* C2PAS: Exit *) Result := y1*(-sin(f1)*cos(f2)+cos(f1)*sin(f2))+y2*(-sin(f1)*sin(f2)-cos(f1)*cos(f2));
  end;
3:
  begin
   (* C2PAS: Exit *) Result := y1*(-cos(f1)*sin(f2)+sin(f1)*cos(f2))+y2*(cos(f1)*cos(f2)+sin(f1)*sin(f2));
  end;
4:
  begin
   (* C2PAS: Exit *) Result := y1*(-cos(f1)*cos(f2)-sin(f1)*sin(f2))+y2*(-cos(f1)*sin(f2)+sin(f1)*cos(f2));
  end;
5:
  begin
   (* C2PAS: Exit *) Result := y1*(sin(f1)*sin(f2)+cos(f1)*cos(f2))+y2*(-sin(f1)*cos(f2)+cos(f1)*sin(f2));
  end;
end;
end;

end.











