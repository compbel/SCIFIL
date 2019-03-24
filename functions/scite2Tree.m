function AM = scite2Tree(addr,n,m)
fid = fopen(addr,'r');
AM = zeros(n+m+1,n+m+1);

line = fgets(fid);
line = fgets(fid);
for i = 1:m
    line = fgets(fid); 
    data = sscanf(line,'%i %s %i;');
    u = data(1);
    v = data(4);
    AM(u,v) = 1;
end
line = fgets(fid);
while true
    line = fgets(fid);
    if strcmp(line(1),'}')
        break;
    end
    data = textscan(line,'%u %s s%u;');
    u = data{1};
    v = data{3}+1;
    AM(u,m+v+1) = 1;
end
ind = [m+1 1:m (m+2):(m+n+1)];
AM = AM(ind,ind);
fclose(fid);
