function [y,id,firmid,id_old,firmid_old,firmid_original,controls,year] = connected_set(y,id,firmid,lagfirmid,firmid_original,controls,year,akm_model_f)
%Finding the largest connected set, using the code from CHK(2013)


%Save ids
firmid_old=firmid;
id_old=id;

%RENAME
N=length(y);
sel=~isnan(lagfirmid);

%relabel the firms
[firms,m,n]=unique([firmid;lagfirmid(sel)]);

firmid=n(1:N);
lagfirmid(sel)=n(N+1:end);


%relabel the workers
[ids,m,n]=unique(id);
id=n;

%initial descriptive stats 
if 0 == 1
s=['Original Dataset:'];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['# of p-y obs: ' int2str(length(y))];
disp(s);
s=['# of workers: ' int2str(max(id))];
disp(s);
if akm_model_f
    s=['# of firms: ' int2str(max(firmid))];
elseif akm_model_fy
    s=['# of firm-years: ' int2str(max(firmid))];
end
disp(s);
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
%FIND CONNECTED SET
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Finding connected set...'];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
end
A=sparse(lagfirmid(sel),firmid(sel),1); %adjacency matrix
%make it square
[m,n]=size(A);
if m>n
    A=[A,zeros(m,m-n)];
end
if m<n
    A=[A;zeros(n-m,n)];
end
A=max(A,A'); %connections are undirected
%[sindex, sz]=components(A); %old code using BGL.
[sindex,sz] = conncomp(graph(A));
sz=sz';
sindex=sindex';
idx=find(sz==max(sz)); %find largest set
%s=['# of firms: ' int2str(length(A))];
%disp(s);
%s=['# connected sets:' int2str(length(sz))];
%disp(s);
%s=['Largest connected set contains ' int2str(max(sz)) ' firms'];
%disp(s);

firmlst=find(sindex==idx); %firms in connected set
sel=ismember(firmid,firmlst);

y=y(sel);
id=id(sel);
firmid=firmid(sel);
firmid_original = firmid_original(sel);
firmid_old=firmid_old(sel);
id_old=id_old(sel);
controls=controls(sel,:);
year=year(sel);

%relabel the firms
[firms,m,n]=unique(firmid);
firmid=n;

%relabel the workers
[ids,m,n]=unique(id);
id=n;
end

