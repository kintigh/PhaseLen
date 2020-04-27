program PhaseLen75;
{$Q+}
{$APPTYPE CONSOLE}
{%File 'ModelSupport\default.txvpck'}

uses
  SysUtils,
  KWKStdXE;
{$R+}
{note this only works for the range of dates in Const from ideal_origin}
{to ideal_limit within the overall limits of dates of the range of integers}

const   version='7.5';  {fix Bugs}
                {enable use not doing the Monte Carlo Analysis}
                {Make sampling of SD's w/o replacement, calculate K properly}
                {7.2 adjust KS output, add CSV input EOF bug}
                {7.1 weighted mean of intercepts for calibrated interval, also Ihat}
                {7.0 increase KS, trials limits added 1ntcal13.14c 2020)
                {5.1  Windows, debugging, IQR}
                {5.0 Calibrated Dates}
        years='1996-2020';
        { to change any of theses values recompile the program }
        ideal_origin=-50000;{-10000}{earliest date considered by the KS method}
        ideal_limit=2500;   {3000}  {latest date considered by the KS method}
        calend=1;            {subscript for calib}
        c14=2;
        c14sd=3;                {these values are taken in but not used}
        calib_year_array_earliest=-50000;  {lower bound of calib_year considered}
        maxidealtrial=100000; {to create ideal distribution for model}
        maxtrial=1000000;  {20000}    {number of date samples drawn}
        maxinterval=200;   {number of phase intervals considered}
        maxdate=2500;      {must accommodate # of empiciacl or test sample
                            dates and number of intercepts of the dates}
        maxtest=5000;      {maximum number of test runs in test mode}

        max_calib_dataset_lines=6000; {works for intCal13}
        rounddates=5;         {round all output dates to the nearest...years}

type    {sample_array=array[1..maxtrial] of real;}  {1..maxtrial}
        date_array=array[1..maxdate] of integer;
        int_array_type=array[1..1] of integer;
        int_array_pntr=^int_array_type;
        pathtype=ansistring;

var     date, date_std, date_intercept, date_index, date_count, null_array, calendar_date:
          date_array;
        ideal: int_array_pntr;

        rslt_span,rslt_len_c,rslt_iqr_c {,rslt_maxdif_c}: array[0..maxinterval] of integer;
        test_array: array [0..maxtest,1..6] of integer;

        rslt_maxdif, rslt_maxdif_mean, rslt_dst, rslt_dst_s,
          rslt_len, rslt_len_p, rslt_len_s, rslt_iqr, rslt_iqr_p, rslt_iqr_s,
          rslt_maxdif_s{, rslt_maxdif_p}: array[0..maxinterval] of real;

        calib: array[0..max_calib_dataset_lines,calend..c14SD] of integer;
        calib_year: int_array_pntr;

        firstdate, lastdate, daterange, mean_date, ndate, median_date,
          {maxempiricalstd,} date_iqr, cdate_iqr, orig_randseed, inc_width,
          cfirstdate,clastdate,cdaterange,cmean_date,cmedian_date, nintercept,
          calib_entries, calib_earliest, calib_latest, calib_origin,
          ninterval,nidealtrials, ntrials, maxspan, middle, usermiddle,
          minspan, testspan, ntrialadj, ideal_date,
          generate_start, generate_finish, generate_length,generate_middle,generate_std,
          sd_loop_start, sd_loop_inc, sd_loop_stop, this_loop, all_loops,
          ndate_loop_start, ndate_loop_inc, ndate_loop_stop,
          e_ihat, test_no, n_tests, ihat_warning: integer;

        start_time,time_now,last_time,elapsed_time,lbound,ubound,{sample_date,}
          ihatk: real;

        calibin, filename, filein:  ansistring;
        dskin,dskout:text;
        median,mode,model: ansichar;
        warn_neg, euclid, big_loop, calibrated_dates, MonteCarlo: boolean;

{************************************}
{**   General purpose utilities    **}
{************************************}

function max(v1,v2: integer): integer;
begin
  if v1<v2 then max:=v2 else max:=v1;
end;

function min(v1,v2: integer): integer;
begin
  if v1<v2 then min:=v1 else min:=v2;
end;

function roundto(v,r: integer): integer;
{Round i to the nearest r}
var i: integer;
begin
  i:=(v mod r); {remainder}
  if i>=((r+1) div 2) then v:=v+(r-i) else v:=v-i;
  roundto:=v;
end;

{Get a truncated normal value between low and up expressed as SD}
function tnormal(low,up: real): real;
var val: real;
begin
  val:=-1.0e99;
  while (val<low) or (val>up) do val:=normal;
  tnormal:=val;
end;

{Sorts a Vector of Dates in Order}
procedure sortdate(n: integer; var vector, col2: date_array);
var h,i,j,s,k,k1,p: integer;  stop: boolean;
begin  { quicksort see knuth Vol 3 }
 {Calculate the quicksort partition size}
  p:=1;
  while (twotothe(p+1)-1)<(n div 3) do inc(p);
  for s:=p downto 1 do begin
    h:=twotothe(s-1);
    for j:=h+1 to n do begin
      stop:=false;  i:=j-h;  k:=vector[j]; k1:=col2[j];
      while (i>0) and not stop do
        if k<vector[i] then begin
          vector[i+h]:=vector[i];
          col2[i+h]:=col2[i];
          i:=i-h;
        end
        else stop:=true;
      vector[i+h]:=k;
      col2[i+h]:=k1;
    end;
  end;
end;

procedure randomize_order(n: integer; var values: date_array);
var i: integer;
    vector, indx, tvalues: date_array;
begin
  for i:=1 to n do begin
    indx[i]:=i;
    vector[i]:=random(100000);
    tvalues[i]:=values[i];
  end;
  sortdate(n,vector,indx);
  for i:=1 to n do values[i]:=tvalues[indx[i]];
end;

{read a set of input dates for evaluation}
procedure readdates;
var cnt,i,j,repetition,nval: integer;
    keyentry: boolean;
begin
  firstdate:=maxint;  lastdate:=0-maxint;  warn_neg:=false;
  filein:='.ADF';
  readfile('File with Empirical Dates (Con for Keyboard)',dskin,filein);
  keyentry:=(filein='CON');
  if not keyentry then begin
    ndate:=getinteger(dskin);
    nval:=getinteger(dskin);

    if (nval<2) or (nval>3) then begin
      writeln('Input Data must have 2 {or 3} values per line, mean, std, {count}');
      halt;
    end;
    cnt:=0;
    for i:=1 to ndate do begin
      inc(cnt);
      date[cnt]:=round(getreal(dskin));
      if date[cnt]<0 then warn_neg:=true;
      date_std[cnt]:=round(getreal(dskin));
      if nval=3 then begin {third value is a count}
        repetition:=getinteger(dskin);
        for j:=2 to repetition do begin
          inc(cnt); date[cnt]:=date[cnt-1];
          date_std[cnt]:=date_std[cnt-1];
        end;
      end;
    end;
    ndate:=cnt;
    close(dskin);
  end
  else {if keyentry then} begin
    ndate:=readint('Number of Actual Dates',1,maxdate,'');
    for i:=1 to ndate do begin
      date[i]:=round(readreal('  Date Mean',-15000,15000,''));
      if date[i]<0 then warn_neg:=true;
      date_std[i]:=round(readreal('  Date Std ',0,10000,''));
    end;
  end;
  sortdate(ndate,date,date_std);
end;

function midspread(n: integer; var datelist: date_array): integer;
{find the midspread of a sorted list}
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

{*************************************}
{****   Calibration Procedures   *****)
{*************************************}

function offset(v,origin: integer): integer;
begin                  {conceptual array with lower bound of origin}
  offset:=v-origin+1;  {offset in actual linear array with bounds from 1..?}
end;

function interpolate(pos,fromss,toss,finddate: integer): integer;
var dif,frac: real;
begin
  dif:=calib[pos+1,fromss]-calib[pos,fromss];
  if dif=0 then writeln('Possible problem in Interpolate, dif=0, year=',finddate);
  if dif=0 then frac:=0
  else frac:=1.0*(finddate-calib[pos,fromss])/dif;
  interpolate:=round(calib[pos,toss]+frac*(calib[pos+1,toss]-calib[pos,toss]));
end;

procedure fill_calib_year;
var i,j,k: integer; {c: char;}
begin
  i:=1;
  if calib_earliest<calib_year_array_earliest then begin
    calib_origin:=calib_year_array_earliest;
   {skip over entries before earliest considered}
    while (i<calib_entries) and (calib_origin>calib[i,calend]) do inc(i);
    calib_origin:=calib[i,calend];
  end
  else calib_origin:=calib_earliest;

  alloc_matrix(calib_year,calib_latest-calib_origin+1,sizeof(integer));

  {fill calib_year with corresponding interpolated c14 dates
   calib_year index is a BP calendar year; value is the
   linearly interpolated radiocarbon year between adjacent values
   in calibration dataset}

  {$R-}
  calib_year^[offset(calib[i,calend],calib_origin)]:=calib[i,c14];
  for j:=i to calib_entries-1 do
    for k:=calib[j,calend]+1 to calib[j+1,calend] do
      calib_year^[offset(k,calib_origin)]:=interpolate(j,calend,c14,k);
  {for i:=calib_origin to calib_latest do begin
    write(i:6,offset(i,calib_origin):7,calib_year^[offset(i,calib_origin)]:7);
    if (i mod 88)=0 then begin
      c:=readchoice('[C]ontinue [Q]uit [S]kip?','CQ','C');
      if c='Q' then halt;
    end;
  end;}
  {$R+}
end;

procedure read_calib;
var c,df: ansichar; i: integer; lne: ansistring;
    calendar_year, cal_BP, skip, radiocarbon, sigma: real;
begin
  {if date[ndate]>6000 then df:='2' else df:='1';} df:='1';
  c:=readchoice('Calibration File: [1]IntCal13, [2]IntCal93, [3]UWTen93, [M]arine93','123M',df);
  case c of
  '1': calibin:='INTCAL13.14C';
  '2': calibin:='INTCAL93.14C';
  '3': calibin:='UWTEN93.14C';
  'M': calibin:='MARINE93.14C';
  end;
  readfile('Calibration File',dskin,calibin);
  {writeln('  Note Program will only work on calibrated dates back to ',
    0-calib_year_array_earliest,' BP');}         {not clear why this limit}
   for i:=1 to 11 do begin
    readln(dskin,lne); {skip over header}
    {writeln(lne);}       {for debugging}
  end;
  calib_entries:=0;
  calib[0,1]:=low(integer); calib[0,2]:=maxint; calib[0,3]:=maxint;

  while not(seekeof(dskin)) do begin
    inc(calib_entries);
    case c of
     '1': begin   {IntCal13.14c}
        cal_BP:=getreal(dskin);
        radiocarbon:=getreal(dskin);
        skip:=getreal(dskin);    {error}
        skip:=getreal(dskin);    {delta14C}
        sigma:=getreal(dskin);   {sigma}
        if cal_BP<1950 then calendar_year:=1950-cal_BP else calendar_year:=1950-cal_BP-1;
      end;
      '2','3','M': begin
        {readln(dskin,calendar_year,delta14C,sigma,radiocarbon,error);}
        calendar_year:=getreal(dskin);
        skip:=getreal(dskin);   {delta14C}
        sigma:=getreal(dskin);       {sigma}
        radiocarbon:=getreal(dskin);
        skip:=getreal(dskin);
      end;
    end;

    if calendar_year>10000 {gets rid of 99999.9 for missing data} then dec(calib_entries)
    else begin
      calib[calib_entries,calend]:=round(calendar_year);
      calib[calib_entries,c14]:=round(radiocarbon);
      calib[calib_entries,c14SD]:=round(sigma);
    end;
  end;
  close(dskin);
  calib_earliest:=calib[1,calend];
  calib_latest:=calib[calib_entries,calend];
  fill_calib_year;
end;

function next_intercept(var pos: integer; finddate: integer): integer;
var rslt: integer;
begin
  while (pos<calib_entries) and (finddate<>calib[pos,c14]) and
    (not ((finddate<max(calib[pos,c14],calib[pos+1,c14])) and
          (finddate>min(calib[pos,c14],calib[pos+1,c14])))) do inc(pos);
  if pos=calib_entries then rslt:=maxint else
  begin
    if finddate=calib[pos,c14] then rslt:=calib[pos,calend]
    else rslt:=interpolate(pos,c14,calend,finddate);
    inc(pos);
  end;
  next_intercept:=rslt;
end;

procedure calendar_stats(count: integer; var date_vector,weight: date_array; weighted: boolean);
var weightsum, i: integer; sum: real;
begin
  sortdate(count,date_vector,weight);

  cfirstdate:=date_vector[1];
  clastdate:=date_vector[count];
  cdaterange:=clastdate-cfirstdate;
  cdate_iqr:=midspread(count,date_vector);

  if (count mod 2)=1 then cmedian_date:=date_vector[(count+1) div 2] else
  cmedian_date:=round((1.0*date_vector[(count div 2)]+1.0*date_vector[(count div 2)+1])/2);

  sum:=0.0; weightsum:=0;
  for i:=1 to count do
    if weighted then begin
      sum:=sum+date_vector[i]/weight[i];
      weightsum:=weightsum+weight[i];
    end
    else begin
      sum:=sum+date_vector[i];
      inc(weightsum);
    end;
  cmean_date:=roundto(round(sum/weightsum),rounddates);

end;

procedure calibrate;
{take uncalibrated dates in date and put calibrated intercepts in datemean}
var i,j,d, pos,start{,cnt}: integer;  {sum: real; }
begin
  nintercept:=1;
  for i:=1 to ndate do begin
    pos:=1;
    while pos<>calib_entries do begin
      date_intercept[nintercept]:=next_intercept(pos,date[i]);
      date_index[nintercept]:=i;
      if date_intercept[nintercept]<>maxint then inc(nintercept);
      if nintercept>maxdate then begin
        HaltProgram('Fatal Error: Problem too Large: Too Many Calibrated Intercepts');
      end;
    end;
  end;
  dec(nintercept);
  if nintercept<=0 then begin
    HaltProgram('Fatal Error: All Dates out of range of calibration file');
  end;
  {this is put in to weight each real date equally and weight the intercepts accordingly}
  i:=1;
  d:=1;
  start:=1;

  while i<=nintercept do begin
    while (i<=nintercept) and (date_index[i]=d) do inc(i);
    for j:=start to (i-1) do date_count[j]:=i-start;
    inc(d); start:=i;
  end;
end;

{*************************************}
{****    Date Calculations       *****)
{*************************************}

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

function calculate_ihat: integer;
var i,ihat: integer;
    meandate, meanstd, datesum, stdsum, datesumsq, datevar, adjvar, avgvar: real;
begin
  if ndate<=1 then ihat:=-1
  else begin
    datesum:=0.0; datesumsq:=0.0;   stdsum:=0.0;
    for i:=1 to ndate do begin
      datesum:=datesum+date[i];
      datesumsq:=datesumsq+1.0*date[i]*date[i];
      stdsum:=stdsum+date_std[i];
    end;
    meandate:=datesum/ndate;
    meanstd:=stdsum/ndate;
    datevar:=(1.0/(ndate))*(datesumsq-1.0*ndate*meandate*meandate) ;

    adjvar:=datevar*ndate/(ndate-1);
    avgvar:=1.0*meanstd*meanstd;

    if adjvar>=avgvar then begin
      ihat:=round(ihatK*sqrt(adjvar-avgvar));
    end
    else begin
      ihat:=-1;
      inc(ihat_warning);
    end
  end;
  calculate_ihat:=ihat;
end;

procedure date_stats;
var i: integer; total: real;
begin;
  {maxempiricalstd:=0; }
  total:=0.0;
  for i:=1 to ndate do begin
    date_intercept[i]:=date[i];  {Sorted vector of original dates}
    total:=total+date[i];
    {if date_std[i]>maxempiricalstd then maxempiricalstd:=date_std[i];}
  end;

  firstdate:=date[1];
  lastdate:=date[ndate];
  daterange:=lastdate-firstdate;
  date_iqr:=midspread(ndate,date);

  mean_date:=round(total/ndate);
  if (ndate mod 2)=1 then median_date:=date[(ndate+1) div 2] else
  median_date:=round((1.0*date[(ndate div 2)]+1.0*date[(ndate div 2)+1])/2);
  e_Ihat:=calculate_Ihat;
end;

{*************************************}
{****   Input/Output Procedures  *****)
{*************************************}

procedure write_parms_console(head: boolean);
begin
  if head then begin
    writeln;
    writeln({'Random Seed ',}'#Dates  Earliest  Latest   Span    IQR   Mean  Median  Ihat    K');
  end;
  write({orig_randseed:11,' ',}ndate:6,firstdate:10,lastdate:8,
    daterange:7,date_iqr:7,mean_date:7,median_date:8,e_ihat:6,ihatk:5:2,' Uncalibrated');
  if e_ihat<0 then writeln(' (Ihat undefined)') else writeln;
  if calibrated_dates then
    writeln({' ':12,}nintercept:6,cfirstdate:10,clastdate:8,
      cdaterange:7,cdate_iqr:7,cmean_date:7,cmedian_date:8,'            ','BC/AD');
end;

procedure print_date_stats(head:boolean);
begin
  if head then begin
    writeln(dskout);
    writeln(dskout,'#Dates  Earliest  Latest   Span    IQR   Mean  Median  Ihat    K');
  end;
  write(dskout,ndate:6,firstdate:10,lastdate:8,
    daterange:7,date_iqr:7,mean_date:7,median_date:8,e_ihat:6,ihatk:5:2,' Uncalibrated');
  if e_ihat<0 then writeln(dskout,' (Ihat is undefined)')else writeln(dskout);
  if calibrated_dates then
    writeln(dskout,nintercept:6,cfirstdate:10,clastdate:8,
      cdaterange:7,cdate_iqr:7,cmean_date:7,cmedian_date:8,' ':12,'BC/AD');
  {writeln(dskout);}
end;

procedure output_date_stats(headonly:boolean);
begin
  if headonly then begin
    write(dskout,'Dates, True_Start, True End, Model, Model_SD, Earliest, Latest, Span, IQR, Mean, Median, Ihat, K');
    if calibrated_dates then
      writeln(dskout,', Intercepts, Cal_Earliest, Cal_Latest, Cal_Span, cal_IQR, Cal_Mean, Cal_Median')
    else writeln(dskout)
  end else begin
    write(dskout,ndate,',',generate_start,',',generate_finish,', "',model,'",',ubound:5:2,',',
      firstdate,',',lastdate,',',daterange,',',date_iqr,',',mean_date,',',
      median_date,',',e_ihat:6,',',ihatk:5:2);
    if calibrated_dates then
      writeln(dskout,',',nintercept,',',cfirstdate,',',clastdate,',',
      cdaterange,',',cdate_iqr,',',cmean_date,',',cmedian_date)
    else writeln(dskout);
  end;
end;

procedure write_test_results(var fle: text; extension: ansichar);
begin
  if MonteCarlo then
    Case extension of
      'S': writeln(fle,n_tests*(this_loop-1)+test_no:4,
        this_loop:4,test_no:5,test_array[0,1]:8,rslt_span[test_array[0,2]]:8,
        rslt_span[test_array[0,3]]:8,rslt_span[test_array[0,4]]:8,
        rslt_span[test_array[0,5]]:8,test_array[0,6]:8,
        ndate:5,generate_std:5,' "N',ndate:3,'/S',generate_std:3,'"');
      'C': writeln(fle,n_tests*(this_loop-1)+test_no,',',
        this_loop,',',test_no,',',test_array[0,1],',',rslt_span[test_array[0,2]],',',
        rslt_span[test_array[0,3]],',',rslt_span[test_array[0,4]],',',
        rslt_span[test_array[0,5]],',',test_array[0,6],',',
        ndate,',',generate_std:5,',',' "N',ndate:3,'/S',generate_std:3,'"');
    end
  else
    Case extension of
      'S': begin
             write(fle,ndate:6,firstdate:10,lastdate:8,
               daterange:7,date_iqr:7,mean_date:7,median_date:8,e_ihat:6,ihatk:5:2);
             if calibrated_dates then
               writeln(fle,nintercept:6,cfirstdate:10,clastdate:8,
               cdaterange:7,cdate_iqr:7,cmean_date:7,cmedian_date:8)
             else writeln(fle);
           end;
      'C': begin
             write(fle,n_tests*(this_loop-1)+test_no,',',
               this_loop,',',test_no,',',firstdate,',',lastdate,',',
               daterange,',',date_iqr,',',mean_date,',',median_date,',',e_ihat,',',ihatk:5:2,',',
               ndate:5,',',generate_std:5,',',' "N',ndate:3,'/S',generate_std:3,'"');
             if calibrated_dates then
               writeln(fle,',',nintercept,',',cfirstdate,',',clastdate,',',
               cdaterange,',',cdate_iqr,',',cmean_date,',',cmedian_date)
             else writeln(fle);
           end;
    end;
end;

procedure print_run_parms;
begin
  Writeln(dskout,'Program PhaseLen V',version,' - (c) ',years,
    ' Keith W. Kintigh');
  writeln(dskout);
  writeln(dskout,'File: ',filename);
  writeln(dskout,'Random Number Seed:       ',orig_randseed);
  if mode<>'T' then print_date_stats(true);

  if montecarlo then begin
    writeln(dskout,'KS Population Size:       ',nidealtrials);
    writeln(dskout,'Sample Comparison Trials: ',ntrials);
  end;

  if calibrated_dates then writeln(dskout,'Calibrated (Calendar Year) Interval Used')
  else writeln(dskout,'Uncalibrated Interval Used');

  if model='R' then writeln(dskout,'Model Distribution: Rectangular')
  else writeln(dskout,'Model Distribution: Truncated Normal +/- ',
    ubound:5:2,' sd');

  if montecarlo then begin
    if median='2' then writeln(dskout,'Using Median as Center of True Intervals')
     else if median='1' then writeln(dskout,'Using Mean as Center of True Intervals')
        else writeln(dskout,' Using User-specified Midpoint of ',usermiddle:5,
          ' as Center of True Intervals');
    if euclid then writeln(dskout,'Distance = Euclidean Distance/Number of Dates')
      else writeln(dskout,'Distance= Mean|deviation|');
  end;
end;

procedure model_dist_parms;
begin
  writeln('Model Distribution of True Dates Across the Interval');
  model:=readchoice('  Model: [R]ectangular or Truncated [N]ormal','RN','R');
  if model='N' then begin
    writeln('  Increasing the std. cutoff increases the weight of the center');
    {e.g. 1 makes -1=phase begin, +1=phase end}
    lbound:=-readreal('  Std. Cutoff: +/-',0.1,5.0,'2.0');
    ubound:=-lbound;
  end else model:='R';
  ihatK:=calc_K(model='R',ubound);
  {Asymmetric unimplemented because it is no longer centered}
  {else if model='A' then  begin
    writeln('  Asymmetric note: 0,1 gives decline, -1, 0 increase over phase');
    lbound:=readreal('  Low Std Cutoff',-10.0,10.0,'-1.0');
    ubound:=readreal('  Up Std Cutoff ',lbound+0.1,10.0,'1.0');
  end};
end;

procedure readparm;
var  dflt,temp: string[6]; datespan,t:integer; {k: real;}
begin
  filename:=file_prefix(filein)+'.TXT';
  writefile('Listing File or Device',dskout,filename);

  if mode<>'T' then model_dist_parms;
  {e_Ihat:=calculate_Ihat; }

  if montecarlo then
  begin
    {nidealtrials:=100000; }
    if mode='T' then dflt:='10000' else dflt:='10000';
    nidealtrials:=readint('Number of Trials to Create KS Population',1,maxidealtrial,dflt);

    if mode='T' then dflt:='1000' else dflt:='10000';
    ntrials:=readint('Number of Trials for Sample Comparisons',0,maxtrial,dflt);
    if ntrials=1 then ntrials:=2;

    median:=readchoice('Interval Midpoint Set at [1] Mean, [2] Median, [U]ser Selected Midpoint','12U','1');
    if median='U' then usermiddle:=readint(' User Specificed Midpoint of True Intervals',-maxint,maxint,'');

    euclid:=readchoice('Use Adjusted [E]uclidean Distance or Mean|[D]eviation|',
      'ED','E')='E';

    if calibrated_dates then datespan:=cdaterange else datespan:=daterange;
    minspan:=readint('Minimum True Interval Considered',0,4*datespan,'0');

    datespan:=roundto(datespan+4,10);   {round up to nearest 10}
    if datespan>(maxint div 2) then datespan:=roundto((maxint div 2)-5,10);
    {to avoid problems with very long spans and 2*datespan}
    str(datespan,temp);
    maxspan:=readint('Maximum True Interval Considered',0,4*datespan,temp);
    testspan:=maxspan-minspan;

    if testspan<= 200 then t:=10
       else if testspan<=500 then t:=25
          else if testspan<=1000 then t:=50
             else t:=100;

    str(t,temp);
    inc_width:=readint('Interval Increment (>=2)',2,1000,temp);
    if (testspan mod inc_width)=0 then ninterval:=testspan div inc_width
    else ninterval:=(testspan div inc_width)+1;
    writeln;
  end;
  {e_ihat:=calculate_Ihat;}   {taken out for testmode, but probably not needed fo others}
  if mode<>'T' then write_parms_console(true);
  if filename<>'CON' then print_run_parms;

end;

procedure print_intervals;
var incr: integer;  c1,c2,c3,c4: char;
begin
  writeln(dskout);
  writeln(dskout,
 '          Phase      dmax      Dist     Span (',daterange:4,')  IQR (',date_iqr:4,')');
  writeln(dskout,
 ' Phase -----------  ------ ------------ -----','------', ' -----------');

  writeln(dskout,
 'Length  From    To  Value    Mean   Std  Mean','  %ile', '  Mean  %ile');

  for incr:=0 to ninterval do begin
    if test_array[0,2]=incr then c1:='<' else c1:=' ';
    if test_array[0,3]=incr then c2:='<' else c2:=' ';
    if test_array[0,4]=incr then c3:='<' else c3:=' ';
    if test_array[0,5]=incr then c4:='<' else c4:=' ';
    writeln(dskout,rslt_span[incr]:6,
      middle-(rslt_span[incr] div 2):6,middle-(rslt_span[incr] div 2)+rslt_span[incr]:6,
      rslt_maxdif[incr]:7:3,c1,{rslt_maxdif_p[incr]:6:1,}
      rslt_dst[incr]:7:2,c2,rslt_dst_s[incr]:5:2,
      rslt_len[incr]:6:0,c3,rslt_len_p[incr]:5:1,
      rslt_iqr[incr]:6:0,c4,rslt_iqr_p[incr]:5:1);
  end;
end;

procedure print_empirical;
{Print Simulated Empirical Distribution for Test Mode}
var i: integer;
begin
  writeln(dskout);
  writeln(dskout,
'============================================================================');
  writeln(dskout);
  writeln(dskout,'Loop: ',this_loop,'  Test: ',test_no,'  Number of Dates: ',ndate,
    '  Std: ',generate_std);
  write(dskout,'True Interval: ',generate_start,' to ',
    generate_start+generate_length,'  Middle: ',generate_middle,'  Length: ',
    generate_length,'  Model: ',model);
  if model='N' then writeln(dskout,' (',ubound:4:1,')') else writeln(dskout);
  writeln(dskout);
  writeln(dskout,'"Empirical" Dates (ndate=',ndate,' std=',date_std[1],'): ');
  for i:=1 to ndate do begin
    write(dskout,date_intercept[i]:6);
    if (i<>ndate) and ((i mod 13)=0) then writeln(dskout);
  end;
  writeln(dskout);
end;

procedure print_test_summary;
var t: integer;
begin
  writeln(dskout);
  writeln(dskout,
'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
  writeln(dskout);
  writeln(dskout,'Test Summary (',this_loop,')  Number of Dates: ',ndate,
    '  Standard Deviation: ',generate_std);
  writeln(dskout,'Ihat could not be calculated ',ihat_warning,' times');
  ihat_warning:=0;
  writeln(dskout);
  writeln(dskout,'Test  Middle  dmax    Dist    Span     IQR    Ihat');
  for t:=1  to n_tests do begin
    writeln(dskout,t:4,test_array[t,1]:8,rslt_span[test_array[t,2]]:8,
      rslt_span[test_array[t,3]]:8,rslt_span[test_array[t,4]]:8,rslt_span[test_array[t,5]]:8,test_array[t,6]:8);
  end;
end;

procedure write_heading(var fle: text; extension: ansichar);
begin  {make output file a systat DATA command file to read the data}
  if montecarlo then
  case extension of
    'S': begin
            if pos('_',filename)>0 then
              writeln('Warning: Illegal Name for Systat File, Remove _ from SYC File Name');
            writeln(fle,'New');
            writeln(fle,'Basic');
            write(fle,'Save ',file_prefix(filename),' /S "True Interval ',
              generate_start,' to ',generate_start+generate_length,';  Middle=');
            if median='2' then write(fle,'Median;  Dist=')
            else if median='1' then write(fle,'Mean;  Dist=')
            else write(fle,'User Midpoint=',usermiddle,';  Dist=');
            if euclid then writeln(fle,'Euclidean/NDate"')
            else writeln(fle,'Mean|deviation|"');
            writeln(fle,'Input Sequence Loop Test Mid_Date dmax Dist Span IQR Ihat Ndate SD Group$');
            writeln(fle,'Run');
          end;
    'C': begin
           writeln(fle,'Sequence, Loop, Test, Mid_Date, dmax, Dist, Span, IQR, Ihat, Ndate, SD, Group');
    end;
  end
  else
  case extension of
    {'S': begin
         end; }
    'C', 'S': Begin
           write(fle,'Sequence, Loop, Test, Earliest , Latest,  Span,IQR,   Mean,  Median, Ihat, K, Ndate, SD, Group');
           if calibrated_dates then writeln(fle,', NIntercept, CalEarliest, CalLatest, Calspan, CalIQR, CalMean, CalMedian')
           else writeln(fle);
         End;
  end;
end;

procedure true_interval_parms;
var temp: integer;  dflt: ansichar;
begin
  dflt:='U';
  calibrated_dates:=readchoice('Generate dates from [C]alibrated Radiocarbon or [U]ncalibrated Interval',
      'CU',dflt)='C';
  if calibrated_dates then begin
    ndate:=1;  date[ndate]:=1000;
    read_calib;
  end;

  generate_start:=readint('True Interval Start Date',-100000,100000,'');
  generate_finish:=readint('True Interval End Date',-100000,100000,'');
  {accounts for entering ages BP or dates}
  if generate_finish<generate_start then begin
    temp:=generate_finish;  generate_finish:=generate_start;  generate_start:=temp;
  end;
  generate_length:=generate_finish-generate_start;
  generate_middle:=round((generate_start+generate_finish)/2);
  model_dist_parms;
end;

procedure date_N_SD_parms;
begin
  ndate:=readint('Number of Dates/Sample',1,maxint,'');
  generate_std:=readint('Standard Error of Dates',0{1},1000,'');
end;

{*************************************}
{****   Monte Carlo Functions    *****)
{*************************************}

{Obtain a random date starting with start for an interval of length with a std amd model of m}
{returns calibrated date, if calibrated also return calendar  date }
function get_date(calib: boolean; start,length,std: integer; var calendar: integer; m: ansichar): integer;
var sample: integer;
begin
  calendar:=maxint;
  sample:=maxint;
  case m of
  'R': begin
         {$R-}
         if calib then begin
           calendar:=start+random(length);
           sample:=calib_year^[calendar-calib_origin+1] {offset removed for time}
         end
         {$R+}
         else sample:=start+random(length); {pop. date uniform w/i range}
         sample:=sample+round(normal*std);
       end;
  'N': begin
         {$R-}
         if calib then begin
           calendar:=start+round((ubound+tnormal(lbound,ubound))/(ubound-lbound)*length);
           sample:=calib_year^[calendar-calib_origin+1] {offset removed for time}
         end
         {$R+}
         else sample:=start+round((ubound+tnormal(lbound,ubound))/(ubound-lbound)*length);
         sample:=sample+round(normal*std)
       end;
  end;
  get_date:=sample;
end;

procedure get_maxdif(pop_span, incr:integer);
{create population distribution fpr interval size & find empricial maxdix}
var i,err,start_date,sampledate, ignore: integer;
  sum, dif, maxdif: real;
begin
  sum:=0.0;  err:=0;  ignore:=0;
  start_date:=middle-(pop_span div 2);  {start date of phase}
  for i:=ideal_origin to ideal_limit do
    {$R-} ideal^[offset(i,ideal_origin)]:=0;  {$R+}
  for i:=1 to nidealtrials do begin
    sampledate:=get_date(calibrated_dates, start_date,pop_span,date_std[random(ndate-1)+1],
      ignore,model);
    sum:=sum+sampledate;
    if (sampledate>=ideal_origin) and (sampledate<=ideal_limit) then
      {$R-} inc(ideal^[sampledate-ideal_origin+1]) {$R+}
    else
     inc(err);
  end;
  if err>0 then writeln('Error: ', err,' Dates out of range');

  ideal_date:=round(sum/nidealtrials);

  ntrialadj:=nidealtrials-err;
  if ntrialadj=0 then begin
    writeln('Unable to do KS Comparison - Ideal Array Out of Range');
  end;

  for i:=ideal_origin+1 to ideal_limit do
    {$R-}
    ideal^[i-ideal_origin+1]:=
      ideal^[(i-1)-ideal_origin+1]+ideal^[i-ideal_origin+1]; {offset removed}
    {$R+}

  maxdif:=0;
  if ntrialadj>0 then
    for i:=1 to ndate do begin
      {$R-}
      dif:=abs(1.0*i/ndate-1.0*ideal^[offset(date_intercept[i],ideal_origin)]/ntrialadj);
      {$R+}
      if dif>maxdif then maxdif:=dif;
    end;
  rslt_maxdif[incr]:=maxdif;
end;

{obtain a sample set of dates}
procedure sample_set(pop_span,middle: integer);
const replacement=false;
var i,start_date, ignore: integer;
begin
  randomize_order(ndate,date_std);
  start_date:=middle-(pop_span div 2);  {start date of phase}
  for i:=1 to ndate do begin
    if replacement {stds assigned with replacement} then
      date[i]:=get_date(calibrated_dates,start_date,pop_span,date_std[random(ndate-1)+1],
       calendar_date[i],model)
    else {stds assigned without replacement}
      date[i]:=get_date(calibrated_dates,start_date,pop_span,date_std[i],calendar_date[i], model)
  end;
  sortdate(ndate,date,date_std);
end;

{evaluate the sorted sample stored in date against empirical or test data}
procedure sample_eval(incr: integer);
var i,first,last,span,iqr: integer;  dst, maxdif,dif: real;
begin
  {for range}
  first:=date[1];
  last:=date[ndate];
  span:=last-first;
  rslt_len[incr]:=rslt_len[incr]+span;
  rslt_len_s[incr]:=rslt_len_s[incr]+sqr(1.0*span);
  if span<=daterange then inc(rslt_len_c[incr]);

  {for midspread}
  iqr:=midspread(ndate,date);
  rslt_iqr[incr]:=rslt_iqr[incr]+iqr;
  rslt_iqr_s[incr]:=rslt_iqr_s[incr]+sqr(iqr);
  if iqr<=date_iqr then inc(rslt_iqr_c[incr]);

  {for maxdif}
  maxdif:=0;
  for i:=1 to ndate do begin
    {$R-}
    dif:=abs(1.0*i/ndate-1.0*ideal^[date[i]-ideal_origin+1]/ntrialadj); {offset}
    {$R+}
    if dif>maxdif then maxdif:=dif;
  end;
  {if rslt_maxdif[incr]>maxdif then inc(rslt_maxdif_c[incr]); } {percentile doesn't help}
  rslt_maxdif_mean[incr]:=rslt_maxdif_mean[incr]+maxdif;
  rslt_maxdif_s[incr]:=rslt_maxdif_s[incr]+sqr(maxdif);

  {for dist}
  dst:=0.0;              {deviation from ordered empirical dates}
  if euclid then begin
    for i:=1 to ndate do dst:=dst+sqr(1.0*date[i]-date_intercept[i]);
    dst:=sqrt(dst)/ndate;
  end
  else begin
    for i:=1 to ndate do dst:=dst+1*abs(date[i]-date_intercept[i]);
    dst:=dst/ndate;
  end;
  rslt_dst[incr]:=rslt_dst[incr]+dst;
  rslt_dst_s[incr]:=rslt_dst_s[incr]+sqr(dst);
end;

{Do evaluation trials for Evaluation or Test Modes from min by inc_w for nint intervals}
procedure runtrials(min,inc_w,nint: integer);
var i,incr,{trial,}trialspan,hours,minutes,seconds: integer; variance, percent,time_left: real;
begin
  if mode<>'T' then begin
    writeln;
    writeln(
  '          Phase      dmax      Dist     Span (',daterange:4,')  IQR (',date_iqr:4,')');
  writeln(
  ' Phase -----------  ------ ------------ -----','------',     ' -----------');
    writeln(
  'Length  From    To  Value   Mean   Std  Mean','  %ile',     '  Mean  %ile');
 end;

 for incr:=0 to nint do begin
    trialspan:=min+incr*inc_w;

    rslt_span[incr]:=trialspan;
    get_maxdif(trialspan,incr);

    {sample_date:=0.0; }
    rslt_maxdif_mean[incr]:=0.0;  rslt_maxdif_s[incr]:=0.0;
    {rslt_maxdif_p[incr]:=0.0;     rslt_maxdif_c[incr]:=0; }
    rslt_dst[incr]:=0.0;          rslt_dst_s[incr]:=0.0;
    rslt_len[incr]:=0.0;          rslt_len_s[incr]:=0.0;
    rslt_len_p[incr]:=0.0;        rslt_len_c[incr]:=0;
    rslt_iqr[incr]:=0.0;          rslt_iqr_s[incr]:=0.0;
    rslt_iqr_p[incr]:=0.0;        rslt_iqr_c[incr]:=0;
    for i:=1 to ntrials do begin
      sample_set(trialspan,middle);
      sample_eval(incr);
    end;
    if ntrials>0 then begin
      {sample_date:=sample_date/(1.0*ntrials*ndate); }

      rslt_maxdif_mean[incr]:=rslt_maxdif_mean[incr]/ntrials;
      {rslt_maxdif_s[incr]:=sqrt((rslt_maxdif_s[incr]-
        1.0*ntrials*sqr(rslt_maxdif_mean[incr]))/(ntrials-1));}
      {rslt_maxdif_p[incr]:=100.0*rslt_maxdif_c[incr]/ntrials;}

      rslt_dst[incr]:=rslt_dst[incr]/ntrials;
      variance:=(1/(ntrials-1)) * (rslt_dst_s[incr]-1.0*ntrials*sqr(rslt_dst[incr]));
      if variance<0 then rslt_dst_s[incr]:=variance
      else rslt_dst_s[incr]:=sqrt(variance);

      rslt_len[incr]:=rslt_len[incr]/ntrials;
      rslt_len_s[incr]:=sqrt((rslt_len_s[incr]-
        1.0*ntrials*sqr(1.0*rslt_len[incr]))/(ntrials-1));
      rslt_len_p[incr]:=100.0*rslt_len_c[incr]/ntrials;

      rslt_iqr[incr]:=rslt_iqr[incr]/ntrials;
      rslt_iqr_s[incr]:=sqrt((rslt_iqr_s[incr]-
        1.0*ntrials*sqr(1.0*rslt_iqr[incr]))/(ntrials-1));
      rslt_iqr_p[incr]:=100.0*rslt_iqr_c[incr]/ntrials;
    end;
    if mode<>'T' then
      writeln(rslt_span[incr]:6,
        middle-(rslt_span[incr] div 2):6,middle-(rslt_span[incr] div 2)+rslt_span[incr]:6,
        rslt_maxdif[incr]:7:3,{rslt_maxdif_mean[incr]:6:3,}
        {rslt_maxdif_s[incr]:6:3,rslt_maxdif_p[incr]:6:1,}
        rslt_dst[incr]:7:2,rslt_dst_s[incr]:6:2,
        rslt_len[incr]:6:0, {rslt_len_s[incr]:6:1,}
        rslt_len_p[incr]:6:1,{ideal_date:6,sample_date:6,}
        rslt_iqr[incr]:6:0, rslt_iqr_p[incr]:6:1)

    else begin
      time_now:=timesec;
      if time_now<last_time then time_now:=time_now+24.0*60*60;
      if (time_now-last_time)>5 then begin
        elapsed_time:=elapsed_time+time_now-last_time;
        percent:=(1.0*(this_loop-1)*n_tests*(nint+1)+
          ((test_no-1)*(nint+1)+incr+1))/(1.0*all_loops*n_tests*(nint+1));
        time_left:=elapsed_time/percent-elapsed_time; {seconds}
        hours:=trunc(time_left/3600.0);
        minutes:=trunc((time_left-3600.0*hours)/60.0);
        seconds:=trunc(time_left-3600.0*hours-60.0*minutes);
        write(this_loop:5,' of ',all_loops:4,test_no:6,':',trialspan:4,'    ',
          100.0*percent:6:1,'%  ÷',hours:3,':',minutes:2,':',seconds:2,'     '^M);
        last_time:=time_now;
      end;
    end;
  end;
end;

{*************************************}
{*****   Generate Procedures   *******}
{*************************************}

procedure generate_dates(output_dates: boolean;
  setno: integer; var ihatn: integer; var ihatsum, ihatsumsq: real);
var i: integer;

begin

  if output_dates then filename:='.ADF' else filename:='.CSV';
  if setno=1 then writefile('File for Output',dskout,filename);
  if output_dates then begin
    write(dskout,ndate,', 2 # Generated Uncalibrated Dates ',generate_start,' to ',
      generate_start+generate_length);
    if calibrated_dates then write(dskout,' Calibrated Interval')
    else  write(dskout,' Uncalibrated Interval');
    writeln(dskout,' Model ',model,' (',ubound:5:2,')',' Randseed: ',orig_randseed);
  end
  else if setno=1 then output_date_stats(true);

  {mean:=0.0; first:=maxint;  last:=-maxint;  sumsq:=0.0;}
  for i:=1 to ndate do begin
    calendar_date[i]:=0;
    date[i]:=get_date(calibrated_dates, generate_start, generate_length, generate_std, calendar_date[i], model);
    date_std[i]:=generate_std;
    if output_dates then writeln(dskout,date[i]:6,date_std[i]:5);
  end;

  sortdate(ndate,date,null_array);
  date_stats;
  if calibrated_dates then begin
    sortdate(ndate,calendar_date,null_array);
    calendar_stats(ndate,calendar_date,null_array,false);
  end;

  if e_ihat>=0 then begin
    inc(ihatn);
    ihatsum:=ihatsum+e_ihat;
    ihatsumsq:=ihatsumsq+e_ihat*e_ihat;
  end;

  if setno=1 then write_parms_console(true)
  else if setno<=50 then write_parms_console(false)
  else if setno=51 then writeln('...');

  if not output_dates then output_date_stats(false);
  if output_dates then close(dskout);
end;

procedure generate_loop;
var ihatn, n_sets,setno: integer;
    ihatsum, ihatsumsq: real;
    dates_not_stats: boolean;
begin
  true_interval_parms;
  date_N_SD_parms;

  nintercept:=ndate;
  dates_not_stats:=readchoice('Produce [S]tatistics for Date Sets or [G]enerate 1 Set of Dates','SG','G')='G';
  if not dates_not_stats then
    n_sets:=readint('Number of Sets of Dates to Generate',1,maxint,'1')
  else n_sets:=1;

  ihatsum:=0.0; ihatsumsq:=0.0;   ihatn:=0;   setno:=1;

  if not dates_not_stats {statsitics only} then
    for setno:=1 to n_sets do
      generate_dates(dates_not_stats,setno,ihatn,ihatsum,ihatsumsq)
  else
  repeat
    generate_dates(dates_not_stats,setno, ihatn,ihatsum,ihatsumsq);
  until not readbool('Generate Another Set With the Same Parameters','F');

  ihatsum:=ihatsum/ihatn;
  if ihatn>1 then begin
    writeln;
    writeln('True Interval Length=',generate_length,'  Ihat mean=',ihatsum:7:2,'  std=',
      sqrt((1.0/(ihatn-1))*(ihatsumsq-ihatn*ihatsum*ihatsum)):6:2,
      '  Date Sets=',ndate,'  Date Sets w/ Ihat Defined=',ihatn);
  end;

  if not dates_not_stats then close(dskout);
end;

{*************************************}
{*****   Evaluate Procedures   *******}
{*************************************}

{find the best fit interval for each similarity measure}
procedure find_best(testno: integer);
var min,i: integer;
begin
  {test_array: middle maxdif RMS dist, ihat}
  {note test_array[i,1-4] are indices, test_array[i,5] is ihat value}
  test_array[testno,1]:=middle;

  min:=0;
  for i:=1 to ninterval do
    if rslt_maxdif[i]<rslt_maxdif[min] then min:=i;
  test_array[testno,2]:=min;

  min:=0;
  for i:=1 to ninterval do
    if rslt_dst[i]<rslt_dst[min] then min:=i;
  test_array[testno,3]:=min;

  min:=0;
  for i:=1 to ninterval do
    if abs(rslt_len[i]-daterange)<abs(rslt_len[min]-daterange) then min:=i;
  test_array[testno,4]:=min;

  min:=0;
  for i:=1 to ninterval do
    if abs(rslt_iqr[i]-date_iqr)<abs(rslt_iqr[min]-date_iqr) then min:=i;
  test_array[testno,5]:=min;
  test_array[testno,6]:=e_ihat;
end;

procedure set_middle;
begin
  if (not calibrated_dates) then begin
    if median='2' then middle:=roundto(median_date,rounddates)
    else if median='1' then middle:=roundto(mean_date,rounddates)
    else middle:=usermiddle;
  end else
    if median='2' then middle:=roundto(cmedian_date,rounddates)
    else if median='1' then middle:=roundto(cmean_date,rounddates)
    else middle:=usermiddle;
end;

procedure evaluate_loop;
var dflt: ansistring;
begin
  test_no:=1; n_tests:=1;
  readdates;
  {if warn_neg then dflt:='U' else dflt:='C';}
  dflt:='U';
  writeln('For radicocarbon years BP input answer C or U; otherwise answer U');
  calibrated_dates:=readchoice('  Estimate [C]alibrated or [U]ncalibrated Intervals',
    'CU',dflt)='C';
  if calibrated_dates then begin
    if warn_neg then
      haltprogram('Warning: Negative radiocarbon ages in input; positive ages BP are expected');
    read_calib; {5}
    calibrate;
    calendar_stats(nintercept,date_intercept,date_count,true);
  end;
  date_stats;
  readparm;

  start_time:=timesec;
  set_middle;
  if MonteCarlo then begin
    runtrials(minspan,inc_width,ninterval);
    find_best(0);
    elapsed_time:=timesec-start_time;
    if filename<>'CON' then print_intervals;
  end;

  writeln;
  writeln('Compute Time: ',elapsed_time/60:7:2,' Minutes');

  if filename<>'CON' then begin
    writeln(dskout);
    writeln(dskout,'Compute Time: ',elapsed_time/60:7:2,' Minutes');
    close(dskout);
  end;
end;

{*************************************}
{*****     Test Procedures     *******}
{*************************************}

procedure test_loop_parms;
begin
  true_interval_parms;
  {if calibrated_dates then read_calib;}
  big_loop:=readchoice(
    'Test for [S]ingle or [M]ultiple Sets of Numbers of Dates & SDs','SM','S')='M';
  if big_loop then begin
    ndate_loop_start:=readint('Starting Number of Dates',1,maxint,'');
    ndate_loop_stop:=readint('  Ending Number of Dates',ndate_loop_start,maxint,'');
    if ndate_loop_start=ndate_loop_stop then ndate_loop_inc:=1 else
    ndate_loop_inc:=readint('  Increment in Number of Dates',1,maxint,'5');
    if ((ndate_loop_stop-ndate_loop_start) mod ndate_loop_inc)<>0 then
      ndate_loop_stop:=ndate_loop_start+ndate_loop_inc*
        (((ndate_loop_stop-ndate_loop_start) div ndate_loop_inc)+1);
    sd_loop_start:=readint('Starting Standard Deviation',1,1000,'');
    sd_loop_stop:=readint('  Ending Deviation',sd_loop_start,1000,'');
    if sd_loop_start=sd_loop_stop then sd_loop_inc:=1 else
     sd_loop_inc:=readint('  Increment in Standard Deviation',1,1000,'10');
    if ((sd_loop_stop-sd_loop_start) mod sd_loop_inc)<>0 then
      sd_loop_stop:=sd_loop_start+sd_loop_inc*
        (((sd_loop_stop-sd_loop_start) div sd_loop_inc)+1);
    n_tests:=readint(
      'Number of Test Runs for Each SD and Number of Dates',
        1,maxtest,'20');
    all_loops:=(((ndate_loop_stop-ndate_loop_start) div ndate_loop_inc)+1)*
      (((sd_loop_stop-sd_loop_start) div sd_loop_inc)+1);
  end
  else begin
    date_N_SD_parms;
    ndate_loop_start:=ndate; ndate_loop_inc:=1; ndate_loop_stop:=ndate;
    sd_loop_start:=generate_std;  sd_loop_inc:=1;  sd_loop_stop:=generate_std;
    n_tests:=readint('Number of Test Runs',1,maxtest,'10');
    all_loops:=1;
  end;
end;

procedure test_loop;
var i: integer;  dskout2: text; fileout: pathtype; outfile: ansichar;  extension: ansistring;
begin
  test_loop_parms;

  ndate:=ndate_loop_stop;  generate_std:=sd_loop_stop;
  {Generate an initial set of dates so that Readparm }
  {has values to work from}
  for i:=1 to ndate do date_std[i]:=generate_std;
  sample_set(generate_length,generate_middle);
  if calibrated_dates then begin
     calendar_stats(ndate,calendar_date,null_array,false);
  end;
  date_stats;

  filein:='';
  readparm;

  if montecarlo then outfile:=readchoice('Output Summary File Type [C]SV or [S]ystat','CS','C') else outfile:='C';
  case outfile of
    'C':  extension:='.CSV';
    'S':  extension:='.SYC';
  end;
  fileout:=file_prefix(filename)+extension;
  writefile('Analysis File for Test Summary',dskout2,fileout);
  write_heading(dskout2,outfile);

  elapsed_time:=0.0;  start_time:=timesec;  last_time:=start_time;
  if montecarlo then begin
    writeln;
    writeln('NDate/SD Loop  Test:Interval  %Done  HHHH:MM:SS Remaining');
  end;
  this_loop:=0;

  while (ndate>=ndate_loop_start) and (ndate>0) do begin
    generate_std:=sd_loop_stop;
    nintercept:=ndate;

    while (generate_std>=sd_loop_start) and (generate_std>0) do
    begin
      inc(this_loop);
      {start_time:=timesec;  last_time:=start_time;}
      for i:=1 to ndate do date_std[i]:=generate_std;
      for test_no:=1 to n_tests do begin
        if not((test_no=1) and (this_loop=1)) then begin
          {first go-around,already done}
          sample_set(generate_length,generate_middle);
          if calibrated_dates then calibrate;
          {pick model "empirical" date set}
          date_stats;
        end;

        print_empirical;
        print_date_stats(true);
        if not montecarlo then
          case outfile of
             'S': ;
             'C': ;
          end;
        set_middle;
        if montecarlo then
        begin
          runtrials(minspan,inc_width,ninterval);
          find_best(0);
          if filename<>'CON' then print_intervals;
          if test_no<=maxtest then for i:=1 to 6 do
            test_array[test_no,i]:=test_array[0,i];
          write_test_results(dskout2,outfile);
        end
        else write_test_results(dskout2,outfile);
      end;
      if montecarlo and (n_tests<maxtest) then print_test_summary;
      {total_elapsed:=total_elapsed+elapsed_time;  elapsed_time:=0;}
      generate_std:=generate_std-sd_loop_inc;
    end;
    ndate:=ndate-ndate_loop_inc;
  end;

  if filename<>'CON' then close(dskout);
  if extension='SYC' then begin
    writeln(dskout2,'~');
    writeln(dskout2,'Exit');
  end;
  if fileout<>'NUL' then close(dskout2);
  writeln;
  if ihat_warning>0 then writeln('Ihat could not be calculated ',ihat_warning,' times');
  writeln;
  writeln('Compute Time: ',elapsed_time/60:5:1,' minutes');
end;

{procedure test_code;
var a, b, c, d: integer;
Begin
  {test randomize order} {
  randomize;
  for b:=1 to 20 do begin
    ndate:=10;
    for a:=1 to ndate do begin
      date_std[a]:=100+a;
    end;
    writeln;
     for a:=1 to ndate do write(date_std[a]:5); writeln;
    randomize_order(ndate,date_std);
    for a:=1 to ndate do write(date_std[a]:5); writeln;
  end;
  closewindow;}


  {test loop to figure time for procedure call overhead}
  {start_time:=timesec;
  for a:=1 to 500000 do begin
    b:=random(1000);  c:=random(1000);  d:=offset(b,c);
  end;
  time_now:=timesec;
  for a:=1 to 500000 do begin
    b:=random(1000);  c:=random(1000);  d:=b-c+1;
  end;
  last_time:=timesec;
  writeln('Ratio calc/proc:',(last_time-time_now)/(time_now-start_time):6:3,
    ' Calc, Proc Time:',last_time-time_now:6:4,time_now-start_time:8:4);
  closewindow;
  haltprogram;
end; }

var a: integer;

begin  { main program block }
  {Initialize}
  ihatk:=-1.0;   ihat_warning:=0;   model:='R';
  for a:=1 to maxdate do null_array[a]:=-1;

  copyright('PhaseLen',version,years,
    'Generate Dates or Estimate Interval Length for Dates w/ Gaussian Errors');
  mode:=readchoice('[E]valuate Dates, [G]enerate Dates, [T]est Mode','EGT','E');
  setrandom;
  orig_randseed:=randseed;

  if mode<>'G' then begin
    MonteCarlo:=readbool('Derive Monte Carlo-Based Interval Length Estimates','T');
    alloc_matrix(ideal,ideal_limit-ideal_origin+1,sizeof(integer));
  end;

case mode of
  'G': generate_loop;
  'E':  evaluate_loop;
  'T': test_loop;
end;

  writeln;
  writeln('Program End');
  CloseWindow;
end.

