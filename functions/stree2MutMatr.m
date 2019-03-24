function M_infer = stree2MutMatr(stree)
% stree: [nextChild haplotype parent label frequency timet fitness]
    cells = (find(cellfun(@isempty, stree(:,2)))).';
    mutations = (find(~cellfun(@isempty, stree(:,2)))).';
    M_infer = zeros(length(cells), length(mutations)-1);
    muts = length(mutations);
    for cell=cells
        m = stree{cell,3};
        while (~isempty(stree{m,3}))
            M_infer(cell-muts, stree{m,4}) = 1;
            m = stree{m,3};
        end
    end
    k=1;
    