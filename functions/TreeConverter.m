classdef TreeConverter < handle
    properties
        btree %field to fill recursively
        index = 1
        % bintree: [leftChild righChild parent label frequency timet fitness]
        % stree: [nextChild haplotype parent label frequency timet fitness]
    end
    methods
        
        function stree = convertToStree(self, tree)
            stree = cell(nIntern+nUnCells+1,7);
            for i = 1:(nIntern+nUnCells+1)
                p = find(AM(:,i) > 0);
                if ~isempty(p)
                    stree{i,3} = p;
                end
                ch = find(AM(i,1:(nIntern+1)) > 0);
                stree{i,1} = ch;
                if i <= nIntern + 1
                    h = find(AM(i,(nIntern+2):end) > 0);
                    if ~isempty(h)
                        stree{i,2} = h;
                    end
                    stree{i,4} = i-1;
                end
            end
        end    
        % converts stree to bintree
        function bintree = convertToBinTreeS(self, stree, obsFreqLeafs)
            newTree = self.subtreeFrequency(stree, obsFreqLeafs, 1 );
            n = length(obsFreqLeafs);
            self.index = 0;
            self.btree(n, 7) = 0;
            self.fill_children_s(newTree, 1, 0);
            bintree = self.btree;
            bintree = bintree(1:self.index,:);
        end
        
        function newTree = subtreeFrequency(self, stree, obsFreqLeafs, node)
            newTree = stree;
            freq = obsFreqLeafs(stree{node, 2});
            for i=stree{node,1}
                newTree = self.subtreeFrequency(newTree, obsFreqLeafs, i);
                if isempty(newTree{node, 5})
                    newTree{node, 5} = newTree{i, 5};
                else
                    newTree{node, 5} = newTree{node, 5} + newTree{i, 5};
                end
            end
            if isempty(newTree{node, 5})
                newTree{node, 5} = freq;
            else
                newTree{node, 5} = newTree{node, 5} + freq;
            end
        end
        
        function fill_children_s(self, stree, node, current_parent)
            if (~isempty(stree{node,1}))
                childrenFreq(1, length(stree{node,1})) = 0;
                for i = 1:length(stree{node,1})
                    childrenFreq(i) = stree{stree{node,1}(i), 5};
                end
                [B, I] = sort(childrenFreq,'descend');
                self.index = self.index + 1;
                % write first children to current index
                parent = self.index;
                self.btree(self.index, 1) = self.index + 1;
                self.btree(self.index, 3) = current_parent;
                self.btree(self.index, 4) = stree{ stree{node,1}(I(1)), 2};
                self.btree(self.index, 5) = B(1);
                self.fill_children_s(stree,  stree{node,1}(I(1)), parent)
                
                self.btree(self.index, 2) = self.index + 1;
                for i=2:length(I)
                    self.fill_children_s(stree,  stree{node,1}(I(i)), parent)
                end
            else
                self.index = self.index + 1;
                self.btree(self.index, 1) = self.index + 1;
                self.btree(self.index, 3) = current_parent;
                self.btree(self.index, 4) = stree{node, 2};
                
                self.index = self.index + 1;
                
            end
            
        end
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
