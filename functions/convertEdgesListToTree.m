function tree = convertEdgesListToTree(addr)
fid = fopen(addr,'r');
line = fgets(fid);
%n <wild_type_freq>
data = sscanf(line,'%i %f');
n = data(1);
freq = data(2);
root = TreeNode;
root.nodes = [];
root.value = 0;
root.nodes = [root];
root.children = [];
tree = root;
root.nodes(n+1) = TreeNode;
root.node_index = 1;

cell = TreeNode;
cell.cell_name = 'wild';
cell.freq = freq;
cell.parent = root;
cell.cell_index = 1;
cell.label = [1];%mock value
root.cells = [root.cells cell];
root.children = [root.children cell];

line = fgets(fid);
while line ~= -1
    data  = sscanf(line, '%i %i %f');
    nn = data(1);
    p = data(2);
    freq = data(3);
    node = root.nodes(nn+1);
    parent = root.nodes(p+1);
    parent.children = [parent.children node];
    node.parent = parent;
    node.value = nn;
    node.node_index = nn+1;
    cell = TreeNode;
    cell.parent = node;
    cell.cell_name = strcat('m',int2str(nn));
    cell.label = [1];% mock value
    cell.freq = freq;
    root.cells = [root.cells cell];
    cell.cell_index = length(root.cells);
    node.children = [node.children cell];
    
    
    line = fgets(fid);
end
fclose(fid);
tree = root;
% root.dispTree

