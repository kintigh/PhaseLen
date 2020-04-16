program calculate_K;
{$Q+}
{$APPTYPE CONSOLE}
{%File 'ModelSupport\default.txvpck'}

uses
  SysUtils,
  KWKStdXE;
{$R+}

const   version='1.0';
        years='2020';

function calc_K(rectangular: boolean; limit: real): real;
{limit is +/- sd oof truncted normal}
var sum,sumsq, slice, x, val, s2, mean, K, factor, length, totwt: real;
    nslice, i,cnt: integer;
begin
  {nslice:=readint('Number of slices',100,10000000,'10000');}
  nslice:=20000;
  if rectangular then k:=sqrt(12.0)
    {length:=readreal('Length tested',1,1000,'100');}
    {slice:=length/nslice; sum:=0.0; sumsq:=0.0;
    for i:=0 to (nslice-1) do begin
      x:=slice*i+slice/2.0;
      sum:=sum+x;
      sumsq:=sumsq+x*x;
    end;
    mean:=sum/nslice;
    s2:= (1.0/(nslice-1))*(sumsq-nslice*mean*mean);
  end}
  else
  begin
    {factor:=readreal('Factor',100,10000000,'10000');}
    factor:=1000000;  {values to 6 decimals w/ 20K slides from 0.5  to 4sd)
    {std of normal distribution between sd cutoffs}
    length:=limit*2;
    slice:=length/nslice; sum:=0.0; sumsq:=0.0; totwt:=0; {area:=0;}
    for i:=0 to (nslice-1) do begin
      x:=slice*i-limit+slice/2.0;
      val:=exp(-x*x/2.0)/sqrt(2*pi);
      cnt:=round(factor*val);
      totwt:=totwt+cnt;
      sum:=sum+x*cnt;
      {area:=area+abs(val*slice); }
      sumsq:=sumsq+x*x*cnt;
    end;
    mean:=sum/totwt;
    s2:= (1.0/(totwt-1))*(sumsq-totwt*mean*mean);    {second term should be 0}
    K:=  sqrt(length*length/s2);
  end;
  calc_K:=K
end;

var model:boolean;  tlimit, K: real;

begin  { main program block }
  copyright('calculate_K',version,years,
    'Examine Ihat Performance');

  repeat
    tlimit:=readreal('Normal Std. Cutoff: +/- (0 for Rectangular Distribution)',0,10.0,'2.0');
    model:= tlimit<=0.0;
    K:=calc_K(model,tlimit);
    writeln('Normal SD Bounds +/-=',tlimit:5:2,'SD     K=',K:9:6);
  until not readbool('Calculate another K','T');

  CloseWindow;
end.

