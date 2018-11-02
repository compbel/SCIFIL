classdef TreeConverter < handle
    properties
        btree %field to fill recursively
        index = 1
        % bintree: [leftChild righChild parent label frequency timet fitness]
    end
    methods
        function bintree = convertToBinTree(self, tree)
            %total number of vertices
            if tree.children(1).value == 'X'
                n = tree.subtree_size(tree) - 1 - tree.subtree_size(tree.children(1)); % children(1) is X - hidden mutations
            else
                n = tree.subtree_size(tree);
            end
            n = n*2;
            self.btree = [];
            self.index = 0;
            self.btree(n, 7) = 0;
            subtreeFreq(self, tree);
            
            childrenFreq(1, length(tree.children)) = 0;
            for i = 1:length(tree.children)
                if isempty(tree.children(i).label) && tree.children(i).value ~= 'X'
                    childrenFreq(i) = tree.children(i).subtree_freq;
                else
                    childrenFreq(i) = -1; %mark cells and move them in the end
                end
            end
            [B, I] = sort(childrenFreq,'descend');
            current_parent = self.index;
            for i = 1:length(I)
                if ~isempty(tree.children(I(i)).label) || tree.children(I(i)).value == 'X'
                    break;
                end
                child_index = self.index + 1;
                fill_children(self, tree.children(I(i)), current_parent); % once we fill all left subtree we go to right subtree
                current_parent = child_index;
                if self.btree(child_index, 1) == 0
                    self.index = self.index + 1;
                    self.btree(child_index, 1) = self.index;
                    self.btree(self.index, 3) = child_index;
                    self.btree(self.index, 5) = tree.children(I(i)).cell_freq;
                else
                    tmp = self.btree(child_index, 1);
                    while self.btree(tmp, 1) + self.btree(tmp, 2) ~= 0
                        tmp = self.btree(tmp,2);
                    end
                    self.btree(tmp, 5) = tree.children(I(i)).cell_freq;
                end
                if i < length(I) && isempty(tree.children(I(i+1)).label) && tree.children(I(i+1)).value ~= 'X'
                    self.btree(child_index, 2) = self.index + 1;
                else
                    self.index = self.index + 1;
                    self.btree(child_index, 2) = self.index;
                    self.btree(self.index, 3) = child_index;
                end
            end
            bintree = self.btree;
            if tree.cell_freq ~= 0
                tmp = 1; % search for for right child
                while (bintree(tmp, 2) ~=0)
                    tmp = bintree(tmp, 2);
                end
                if bintree(tmp, 4) == 0 % if it is a child without mutation
                    bintree(tmp, 5) = tree.cell_freq;
                else
                    % add left empty subtree
                    self.index = self.index + 1;
                    bintree(tmp, 1) = self.index;
                    bintree(self.index, 3) = tmp;
                    self.index = self.index + 1; % add right empty subtree
                    bintree(self.index, 5) = tree.cell_freq;
                    bintree(self.index, 3) = tmp;
                    bintree(tmp, 2) = self.index;
                end
            end
            bintree = bintree(1:self.index,:);
        end
        
        function fill_children(self, node, current_parent)
            self.index = self.index + 1;
            %save current node
            self.btree(self.index, 3) = current_parent;
            self.btree(self.index, 4) = node.value;
            childrenFreq = [];
            if ~isempty(node.children)
                childrenFreq(1, length(node.children)) = 0;
            end
            %find order of children
            
            for i = 1:length(node.children)
                if isempty(node.children(i).label)
                    childrenFreq(i) = node.children(i).subtree_freq;
                else
                    childrenFreq(i) = -1; %mark cells and move them in the end
                end
            end
            [B, I] = sort(childrenFreq,'descend');
            if isempty(B) || B(1) == -1
                return
            end
            current_parent = self.index;
            self.btree(self.index, 1) = self.index + 1;% start with all children nodes
            for i = 1:length(I) % next child will go to right subtree of previous child
                child = node.children(I(i));
                if ~isempty(child.label)
                    break;
                end
                child_index = self.index + 1;
                fill_children(self, child, current_parent);
                current_parent = child_index;
                if self.btree(child_index, 1) == 0
                    self.index = self.index + 1;
                    self.btree(child_index, 1) = self.index;
                    self.btree(self.index, 3) = child_index;
                    self.btree(self.index, 5) = node.children(I(i)).cell_freq;
                else
                    tmp = self.btree(child_index, 1);
                    while self.btree(tmp, 1) + self.btree(tmp, 2) ~= 0
                        tmp = self.btree(tmp,2);
                    end
                    self.btree(tmp, 5) = node.children(I(i)).cell_freq;
                end
                if i < length(I) && isempty(node.children(I(i+1)).label)
                    self.btree(child_index, 2) = self.index + 1;
                else
                    self.index = self.index + 1;
                    self.btree(child_index, 2) = self.index;
                    self.btree(self.index, 3) = child_index;
                end
            end
        end
        
        % fills subtree_freq and cell_freq for all nodes recursively
        function f = subtreeFreq(self, node)
            node.cell_freq = 0;
            node.subtree_freq = 0;
            if node.value == 'X'
                node.subtree_freq = 0;
                f = 0;
                return;
            end
            node.cell_freq = 0;
            f = 0;
            if ~isempty(node.freq)
                f = node.freq;
                return;
            end
            for i = 1:length(node.children)
                if ~isempty(node.children(i).freq)
                    node.cell_freq = node.cell_freq + node.children(i).freq;
                end
                f = f + subtreeFreq(self, node.children(i));
            end
            node.subtree_freq = f;
        end
    end
end
