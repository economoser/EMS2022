function [id,firmid,id_old,firmid_old] = connected_set_Chris(id,firmid)
    % Author: Raffaele Saggio. raffaele.saggio@berkeley.edu
    % EDITED BY CHRIS MOSER
    %
    % ONLY Finding the largest connected set, using the code from
    % CHK(2013). Does nothing else and is disconnected from other KSS
    % procedures.
    %
    % INPUTS:
    % - id: original worker ID
    % - firmid: original firm ID
    %
    % OUTPUTS:
    % - id: newly created worker ID
    % - firmid: newly created firm ID
    % - id_old: original worker ID
    % - firmid_old: original firm ID
    
    %Lagfirmid
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker
    clear gcs
    
    %Save ids
    firmid_old=firmid;
    id_old=id;
    
    %RENAME
%     N=length(y); % OLD (by KSS)
    N=length(id); % NEW (by Chris)
    sel=~isnan(lagfirmid);
    
    %relabel the firms
    [~,~,n]=unique([firmid;lagfirmid(sel)]);

    firmid=n(1:N);
    lagfirmid(sel)=n(N+1:end);
    clear N

    %relabel the workers
    [~,~,n]=unique(id);
    id=n;
    clear n

    %initial descriptive stats 
    if 0 == 1
    s=['Original Dataset:'];
    disp(s)
    s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
    disp(s)
    s=['# of worker-years: ' int2str(length(id_old))];
    disp(s);
    s=['# of workers: ' int2str(max(id))];
    disp(s);
    s=['# of firms: ' int2str(max(firmid))];
    disp(s);
    
    %FIND CONNECTED SET
    s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
    disp(s)
    s=['Finding connected set...'];
    disp(s)
    s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
    disp(s)
    end
    
%     XXX OLD BY KSS:
%     A=sparse(lagfirmid(sel),firmid(sel),1); %adjacency matrix
%     clear lagfirmid sel
%     %make it square
%     [m,n]=size(A);
%     if m>n
%         A=[A,zeros(m,m-n)];
%     end
%     if m<n
%         A=[A;zeros(n-m,n)];
%     end
%     clear m n
    
%     XXX BY CHRIS:
    n_rows = max(lagfirmid(sel));
    n_cols = max(firmid(sel));
    n_max = max(n_rows, n_cols);
    A=sparse(lagfirmid(sel),firmid(sel),1,n_max,n_max); %adjacency matrix
    clear n_rows n_cols n_max
    
    A=max(A,A'); %connections are undirected
    %[sindex, sz]=components(A); %old code using BGL.
    [sindex,sz] = conncomp(graph(A));
    sz=sz';
    sindex=sindex';
    idx=find(sz==max(sz)); %find largest set
    clear sz
    %s=['# of firms: ' int2str(length(A))];
    %disp(s);
    %s=['# connected sets:' int2str(length(sz))];
    %disp(s);
    %s=['Largest connected set contains ' int2str(max(sz)) ' firms'];
    %disp(s);

    firmlst=find(sindex==idx); %firms in connected set
    clear sindex idx
    sel=ismember(firmid,firmlst);
    clear firmlst
    
    firmid=firmid(sel); 
    id=id(sel);
    firmid_old=firmid_old(sel);
    id_old=id_old(sel);
    clear sel

    %relabel the firms
    [~,~,n]=unique(firmid);
    firmid=n;
    clear n

    %relabel the workers
    [~,~,n]=unique(id);
    id=n;
    clear n
end

