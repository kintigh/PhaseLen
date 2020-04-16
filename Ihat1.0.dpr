program Ihat10;
{$Q+}
{$APPTYPE CONSOLE}
{%File 'ModelSupport\default.txvpck'}

uses
  SysUtils,
  KWKStdXE;
{$R+}
{note this only works for the range of dates in Const from ideal_origin}
{to ideal_limit within the overall limits of dates of the range of integers}

const   version='1.0';
        years='2020';

        maxdate=10000;      {must accommodate # of dates and number of intercepts of the dates}
        maxtest=100000;      {maximum number of test runs in test mode}

        max_calib_dataset_lines=6000; {works for intCal13}
        rounddates=5;         {round all output dates to the nearest...years}

type    date_array=array[1..maxdate] of integer;


var     date,date_std, date_mean, date_index, date_count: date_array;
        {rslt_array: sample_array;}

        firstdate, lastdate, daterange, mean_date, ndate, median_date,
        {cfirstdate,clastdate,cdaterange,cmean_date,nintercept, cmedian_date,
        maxempiricalstd,  cdate_iqr,
        inc_width, ninterval,nidealtrials, ntrials, maxspan,}
        date_iqr, datepart,  orig_randseed,
        middle, usermiddle, minspan, testspan, ntrialadj,ideal_date,
        generate_start,generate_std,
        sd_loop_start, sd_loop_inc, sd_loop_stop, this_loop, all_loops,
        ndate_loop_start, ndate_loop_inc, ndate_loop_stop, {calib_entries,
        calib_earliest, calib_latest, calib_origin,} e_ihat: integer;

        true_start,true_end,test_no, n_tests, ihat_warning: integer;

        {start_time,time_now,last_time,elapsed_time,}lbound,ubound,sample_date,
        generate_length,generate_middle,ihatk: real;

        filename:  ansistring;
        dskout:text;
        median,mode,model: ansichar;
        warn_neg, euclid, big_loop, output_on: boolean;

{Sorts a Vector of Dates in Order}
procedure sortdate(n: integer; var vector: date_array; p: integer);
var h,i,j,s: integer;  k: integer;  stop: boolean;
begin  { quicksort see knuth Vol 3 }
  for s:=p downto 1 do begin
    h:=twotothe(s-1);
    for j:=h+1 to n do begin
      stop:=false;  i:=j-h;  k:=vector[j];
      while (i>0) and not stop do
        if k<vector[i] then begin
          vector[i+h]:=vector[i];  i:=i-h;
        end
        else stop:=true;
      vector[i+h]:=k;
    end;
  end;
end;

{Calculate the quicksort partition size}
function partition(n: integer): integer;
var p: integer;
begin
  p:=1;
  while (twotothe(p+1)-1)<(n div 3) do inc(p);
  partition:=p;
end;

function max(v1,v2: integer): integer;
begin
  if v1<v2 then max:=v2 else max:=v1;end;

function min(v1,v2: integer): integer;
begin
  if v1<v2 then min:=v1 else min:=v2;
end;

function roundto(v,r: integer): integer;
var i: integer;
begin
  i:=(v mod r); {remainder}
  if i>=((r+1) div 2) then v:=v+(r-i) else v:=v-i;
  roundto:=v;
end;

{find the midspread of a sorted list}
function midspread(n: integer; var datelist: date_array): integer;
var depth: real;  d: integer;
begin
  depth:=(floor((n+1)/2)+1)/2;
  d:=floor(depth)-1;
  if frac(depth)<0.001 then
    midspread:=datelist[n-d]-datelist[1+d]
  else
    midspread:=
      round(1.0*(datelist[n-d]+datelist[n-d-1]-datelist[1+d]-datelist[2+d])/2.0);
end;



function calculate_ihat: integer; {this gets ihat before it is multiplied by k}
{K isn't available at this point in computation}
var i,ihat: integer;
    meandate, meanstd, datesum, datesumsq, datevar, adjvar, avgvar: real;
begin
  if ndate<=1 then ihat:=-1
  else begin
    datesum:=0.0; datesumsq:=0.0;
    for i:=1 to ndate do begin
      datesum:=datesum+date[i];
      datesumsq:=datesumsq+1.0*date[i]*date[i];
    end;
    meandate:=datesum/ndate;
    meanstd:=generate_std;
    datevar:=(1.0/(ndate-1))*(datesumsq-1.0*ndate*meandate*meandate) ;
    avgvar:=1.0*meanstd*meanstd;

    if datevar>=avgvar then begin
      ihat:=round(ihatK*sqrt(datevar-avgvar));
    end
    else begin
      ihat:=-1;
      inc(ihat_warning);
    end
  end;
  calculate_ihat:=ihat;
end;

function calculate_K(model: ansichar; sdtest: real): real; {s is +/i sd oof truncted normal}
var sum,sumsq, slice, x, val, s2, mean, K, factor, length, totwt {, area}: real;
    nslice, i,cnt: integer;
begin
  {nslice:=readint('Number of slices',100,10000000,'10000');}
  nslice:=20000;
  if model='R' then k:=sqrt(12)
    {length:=readreal('Length tested',1,1000,'100');}
    {slice:=length/nslice; sum:=0.0; sumsq:=0.0;
    for i:=0 to (nslice-1) do begin
      x:=slice*i+slice/2.0;
      sum:=sum+x;
      sumsq:=sumsq+x*x;
    end;
    mean:=sum/nslice;
    s2:= (1.0/(nslice-1))*(sumsq-nslice*mean*mean);
  end} else
  begin
   {factor:=readreal('Factor',100,10000000,'10000');}
   factor:=1000000;  {values to 6 decimals w/ 20K slides from 0.5  to 4sd)
   {std of normal distribution between sd cutoffs}
    {new code}
    length:=sdtest*2;
    slice:=length/nslice; sum:=0.0; sumsq:=0.0; totwt:=0; {area:=0;}
    for i:=0 to (nslice-1) do begin
      x:=slice*i-sdtest+slice/2.0;
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
  calculate_K:=K
end;

procedure date_stats;
var i: integer; total: real;
begin;
  total:=0.0;       {this was not initialed in v5}
  for i:=1 to ndate do begin
    date_mean[i]:=date[i];  {Sorted vector of original dates}
    total:=total+date[i];
  end;

  firstdate:={roundto(}date[1]{,rounddates)};
  lastdate:={roundto(}date[ndate]{,rounddates)};
  daterange:=lastdate-firstdate;
  if output_on then date_iqr:=midspread(ndate,date) else date_iqr:=-1;
  {if (firstdate<ideal_origin) or (lastdate>ideal_limit) then begin
    writeln('Earliest date: ',firstdate,'   Last Date: ',lastdate);
    writeln('This problem is out of the date range for the program.');
    writeln('Dates considered must lie between ',ideal_origin,' and ',ideal_limit,'.');
    haltprogram;
  end;}
  mean_date:={roundto(}round(total/ndate){,rounddates)};
  if (ndate mod 2)=1 then median_date:=date[(ndate+1) div 2] else
  median_date:=round((1.0*date[(ndate div 2)]+1.0*date[(ndate div 2)+1])/2);
  {median_date:=roundto(median_date,rounddates);}
  e_ihat:=calculate_ihat;

end;

{Get a truncated normal value between low and up expressed as SD}
function tnormal(low,up: real): real;
var val: real;
begin
  val:=-1.0e99;
  while (val<low) or (val>up) do val:=normal;
  tnormal:=val;
end;

{Obtain a random date starting with start for an interval of length with a std amd model of m}
function get_date(start,length:real; std: integer; m: ansichar): integer;
var sample: real;
begin
  if m='R' then begin
    sample:=start+random*length; {pop. date uniform w/i range}
    sample:=sample+normal*std;
      {sample date; std picked at random from empricial dates}
  end
  else if m='N' then begin
    sample:=start+(ubound+tnormal(lbound,ubound))/(ubound-lbound)*length;
    sample:=sample+normal*std
    {sample date; std picked at random from empricial dates}
  end;
  get_date:=round(sample);
end;

{obtain a sample set of dates}
procedure sample_set(pop_span,middle: real);
var i: integer; start_date: real;
begin
  start_date:=middle-(pop_span/2);  {start date of phase}
  for i:=1 to ndate do begin
    date[i]:=get_date(start_date,pop_span,generate_std,
      model);
    {sample_date:=sample_date+date[i];}
  end;
  if output_on then sortdate(ndate,date,datepart);
end;

procedure loop_parms;
begin

  true_start:=readint('True Interval Start Date',-100000,100000,'');
  true_end:=readint('True Interval End Date  ',true_start+1,100000,'');
  generate_length:=1.0*(true_end-true_start);
  generate_middle:=(true_start+true_end)/2;

  writeln('Model Distribution of True Dates Across the Interval');
  model:=readchoice('  Model: [R]ectangular or Truncated [N]ormal','RN','R');
  if model='N' then begin
    writeln('  Increasing the std. cutoff increases the weight of the center');
    {e.g. 1 makes -1=phase begin, +1=phase end}
    lbound:=-readreal('  Std. Cutoff: +/-',0.1,10.0,'2.0');
    ubound:=-lbound;
  end else model:='R';
  ihatK:=calculate_K(model,ubound);
  ihat_warning:=0;


  {big_loop:=readchoice(
    'Test for [S]ingle or [M]ultiple Sets of Numbers of Dates & SDs','SM','S')='M';}
  big_loop:=false;
  {if big_loop then begin
    ndate_loop_start:=readint('Starting Number of Dates',1,maxint,'');
    ndate_loop_stop:=readint('  Ending Number of Dates',ndate_loop_start,maxint,'');
    if ndate_loop_start=ndate_loop_stop then ndate_loop_inc:=1 else
    ndate_loop_inc:=readint('  Increment in Number of Dates',1,maxint,'5');
    if ((ndate_loop_stop-ndate_loop_start) mod ndate_loop_inc)<>0 then
      ndate_loop_stop:=ndate_loop_start+ndate_loop_inc*
        (((ndate_loop_stop-ndate_loop_start) div ndate_loop_inc)+1);
    sd_loop_start:=readint('Starting Standard Deviation',1,500,'');
    sd_loop_stop:=readint('  Ending Deviation',sd_loop_start,500,'');
    if sd_loop_start=sd_loop_stop then sd_loop_inc:=1 else
     sd_loop_inc:=readint('  Increment in Standard Deviation',1,1000,'10');
    if ((sd_loop_stop-sd_loop_start) mod sd_loop_inc)<>0 then
      sd_loop_stop:=sd_loop_start+sd_loop_inc*
        (((sd_loop_stop-sd_loop_start) div sd_loop_inc)+1);
    n_tests:=readint(
      'Number of Date Sets for Each SD and Number of Dates',
        1,maxtest,'1000');
    all_loops:=(((ndate_loop_stop-ndate_loop_start) div ndate_loop_inc)+1)*
      (((sd_loop_stop-sd_loop_start) div sd_loop_inc)+1);
  end
  else}
  begin
    ndate:=readint('Number of Dates',1,maxdate,'1000');
    generate_std:=readint('Standard Deviation',0{1},1000,'50');
    datepart:=partition(ndate); {setup for sort}
    ndate_loop_start:=ndate; ndate_loop_inc:=1; ndate_loop_stop:=ndate;
    sd_loop_start:=generate_std;  sd_loop_inc:=1;  sd_loop_stop:=generate_std;
    n_tests:=readint('Number of Date Sets',1,maxtest,'20000');
    all_loops:=1;
  end;
end;

var a, b, c, d, i, ihat_n, sample: integer;
    dflt: ansistring;
    tlimit, ihatsum, ihatmean, sample_ihatK: real;

begin  { main program block }
  copyright('Ihat',version,years,
    'Examine Ihat Performance');

 {repeat
    model:=readchoice('Model [R]ectangular or [N]ormal','RN','N');
    if model='R' then tlimit:=0.0 else tlimit:=readreal('Std. Cutoff: +/-',0.01,10.0,'2.0');
    idealK:=calculate_K(model,tlimit);
    writeln('SD Limit=',tlimit:5:2,'   Ideak K=',idealK:9:6);
  until not readbool('Calculate another K','T'); }

  mode:='T';
  Output_on:=readbool('Write output samples to a file','F');
  filename:='.CSV';
  if output_on then writefile('Listing File or Device',dskout,filename);

  setrandom;
  orig_randseed:=randseed;


  if output_on then writeln(dskout,'Sample, NDate, DateStd, Earliest, Latest, Range, IQR, Mean, Median, Ihat, IhatMean, RunID');
  repeat
    loop_parms;
    ndate:=ndate_loop_stop;  datepart:=partition(ndate);  generate_std:=sd_loop_stop;


    this_loop:=0; ihatsum:=0.0;   ihat_n:=0;

    while (ndate>=ndate_loop_start) and (ndate>0) do begin
      generate_std:=sd_loop_stop;

      while (generate_std>=sd_loop_start) and (generate_std>0) do begin
        inc(this_loop);
        {start_time:=timesec;  last_time:=start_time;}
        for test_no:=1 to n_tests do begin
          sample_set(generate_length,generate_middle);
          date_stats;
          if e_ihat>=0 then begin
            ihatsum:=ihatsum+e_ihat;
            inc(ihat_n);
            ihatmean:=ihatsum/ihat_n;
          end;
          sample:= {n_tests*(this_loop-1)+} test_no;

          if output_on then writeln(dskout,sample:4,',',
            ndate:6,',',generate_std:5,',',firstdate:6,',',lastdate:6,',',daterange:6,',',date_IQR:6,',',mean_date:6,',',
            median_date:6,',',e_ihat:6,',',ihatmean:8:2,
            ' "N',ndate:3,'/S',generate_std:3,'"');

        end;
        generate_std:=generate_std-sd_loop_inc;
      end;
      ndate:=ndate-ndate_loop_inc;  datepart:=partition(ndate);
    end;
    sample_ihatk:=1.0*(true_end-true_start)/ihatmean;
    writeln('Date Sets=',sample,'  Numerical IhatK=',ihatk:7:4,'   Ihat Mean=',ihatmean:8:2);
  until not readbool('Run Again','T');
  if output_on then close(dskout);
  writeln;
  writeln('Program End');
  CloseWindow;
end.

