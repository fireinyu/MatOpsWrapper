classdef MatOpsWrapper < handle
    properties
        mat
        ops
        aug = false
    end

    methods (Static)
        
        function obj = build(seq,vars)
            arguments
                seq
                vars = false;
            end
            if (vars == false)
                vars = symvar(seq);
                disp(vars);
            end
            mat = zeros(size(seq,1)*size(seq,2),length(vars));
            rowNum = 0;
            for i = 1:size(seq,1)
                for j = 1:size(seq,2)
                    rowNum = rowNum + 1;
                    expr = seq(i,j);
                    for varIdx = 1:length(vars)
                        var = vars(varIdx);
                        c = coeffs(expr,var);
                        if (length(c)>=2)
                            mat(rowNum,varIdx) = c(2);
                        end
                    end
                end
            end
            obj = qwe(mat);
        end
        function obj = vander(degree,xArr)
            mat = fliplr(vander(xArr));
            mat = mat(:,1:degree+1);
            obj = qwe(mat);
        end

    end
    
    methods
        function obj = MatOpsWrapper(mat)
            obj.mat = sym(mat);
            obj.ops = {};
      
        end
        function pivArr = pivots(obj)
            copyObj = qwe(obj.mat)
            copyObj.rref();
            m = copyObj.mat;
            numRows = size(m,1);
            numCols = size(m,2);
            pivArr = [];
            for rowNum = 1:numRows
                for colNum = 1:numCols
                    if m(rowNum, colNum) ~= 0
                        pivArr = [pivArr colNum];
                        break;
                    end
                end
            end
        end

        function show(obj)
            if obj.aug == false
                disp(obj.mat)
            else
                separator = NaN(size(obj.mat, 1), 1);
                disp([obj.mat,separator,obj.aug.mat])
            end
        end
        function steps(obj)
            for index = 1:length(obj.ops)
                opr = obj.ops{index};
                disp(sprintf("%d|%d:: %s",index-length(obj.ops)-1,index,opr.toString()));
            end
        end
        function oprStr = operate(obj,operation)
            operation.apply(obj);
            if obj.aug == false
                oprStr = obj.mat;
            else
                operation.apply(obj.aug);
                %separator = NaN(size(obj.mat, 1), 1);
                oprStr = obj.mat;
            end
            

        end
        function oprStr = swapRows(obj,r1,r2)
            operation = ElementaryOperation.makeRowSwap(r1,r2);
            oprStr = obj.operate(operation);
        end
        function oprStr = addRow(obj,r1,r2,c)
            operation = ElementaryOperation.makeRowAddition(r1, r2,c);
            oprStr = obj.operate(operation);
        end
        function oprStr = mulRow(obj,r1,c)
            operation = ElementaryOperation.makeRowMultiplication(r1, c);
            oprStr = obj.operate(operation);
        end
        function oprStr = ref(obj)
            operation = CompoundOperation.makeREF(obj);
            oprStr = obj.operate(operation);
            obj.simp();
        end
        function oprStr = rref(obj)
            operation = CompoundOperation.makeRREF(obj);
            oprStr = obj.operate(operation);
            obj.simp();
            
        end
        function opStr = undo(obj,upto)
            for index = 1:upto
                if isempty(obj.ops)
                    break;
                end
                opr = obj.ops{end};
                obj.ops(end) = [];
                opStr = obj.operate(opr.invert());
                obj.ops(end) = [];
                if obj.aug == false
                else
                    obj.aug.ops(end) = [];
                    obj.aug.ops(end) = [];
                end
            end
            obj.simp();
        end
        function augment(obj,mat)
            obj.aug = qwe(mat);
            for index = 1:length(obj.ops)
                obj.ops{index}.apply(obj.aug)
            end
        end
        function copyMat = freezeMat(obj)
            copyMat = obj.mat;
        end
        function copyOps = freezeOps(obj)
            copyOps = obj.ops;
        end
        function copyAug = freezeAug(obj)
            copyAug = obj.aug;
        end
        function copyObj = freeze(obj)
            copyObj = qwe(obj.freezeMat());
            copyObj.ops = obj.freezeOps();
            copyObj.aug = obj.freezeAug();
        end
        function mop = opsMop(obj)
            compoundOps = CompoundOperation(obj.freezeOps());
            mop = compoundOps.toMop(size(obj.mat,1));
        end
        function simp(obj)
            try
                obj.mat = simplify(obj.mat);
            catch exception
             
            end
            if obj.aug ~= false
                obj.aug.simp();
            end
        end
        function delOps(obj, upto)
            if upto > length(obj.ops)
                upto = length(obj.ops);
            end
            obj.ops = obj.ops(upto+1:length(obj.ops));
        end
        function subs(obj,match,replacement)
            obj.mat = subs(obj.mat,match,replacement);
            if obj.aug ~= false
                obj.aug.subs(match,replacement);
            end
        end
        function mop = invert(obj) 
            oriMop = obj.freeze();
            oriMop.delOps(length(oriMop.ops));
            oriMop.augment(eye(size(oriMop.mat,1)));
            oriMop.rref();
            if size(obj.mat,1) == size(obj.mat,2)
                % total inverse
                mop = oriMop.aug;
            elseif size(obj.mat,1) < size(obj.mat,2)
                % right inverse
                mop = qwe([oriMop.aug.mat; zeros([size(obj.mat,2)-size(obj.mat,1),size(obj.mat,1)])]);
            else
                % left inverse
                mop = qwe(oriMop.aug.mat(1:size(obj.mat,2),:));
            end
        end
        function [l, u] = lu(obj)
            mop = obj.freeze();
            mop.delOps(length(mop.ops));
            mop.ref();
            l = mop.opsMop().invert();
            u = mop;
        end
        function mop = subCol (obj, colNum, with)
            mop = obj.freeze();
            mt = mop.mat;
            ag = mop.aug.mat;
            scalars = mt(:,colNum);
            mt(:,colNum) = [];
            ag = ag - scalars * with;
            mop.mat = mt;
            mop.aug.mat = ag;
        end
        function mop = gram(obj)
            %TODO: use in-built qr
            mop = qwe(obj.mat);
            mop.simp();
            m = mop.mat;
            ortho = m(:,1);
            for i = 2:size(m,2)
                pmat = (ortho'*ortho)^(-1)*ortho';
                col = m(:,i);
                x = pmat * col;
                pcol = ortho * x;
                ncol = simplify(col - pcol);
                ortho = [ortho ncol];
                
            end
            [q r] = qr(ortho,"econ");
            mop.mat = simplify(q);

        end
        
        function obj = col(this)
            obj = qwe(colspace(this.mat));
            % pivs = this.pivots();
            % obj = qwe(this.mat(:,pivs));
        end

        function obj = lnull(this)
            obj = qwe(null(this.mat'));
            % rank = length(this.pivots());
            % copyObj = this.freeze();
            % copyObj.au(eye(size(this.mat,1)));
            % copyObj.rr();
            % obj = qwe(transpose(copyObj.aug.mat(rank+1:size(this.mat,1),:)));
        end

        function obj = row(this)
            obj = qwe(colspace(this.mat'));
            % obj = qwe(transpose(this.mat)).col();
        end

        function obj = null(this)
            obj = qwe(null(this.mat));
            % obj = qwe(transpose(this.mat)).lnull();
        end

        function obj = union(this, other)
            obj = qwe([this.mat other.mat]).col();
        end

        function obj = intersect(this, other)
            obj = qwe([this.lnull().mat other.lnull().mat]).lnull();
        end

        function [u, s, v] = svd(this)
            ata = this.mat' * this.mat;
            syms x;
            eigenvalues = sort(solve(det(x*eye(size(ata,1))-ata)),"descend");
            s = sym(zeros(size(this.mat,1),size(this.mat,2)));
            v = sym([]);
            seen = [];
            
            for i = 1:size(eigenvalues,1)
                if (eigenvalues(i)) ~= 0
                    s(i,i) = sqrt(eigenvalues(i));
                end
                if any(seen == eigenvalues(i))
                    continue
                end
                seen = [seen eigenvalues(i)];
                v = [v qwe(ata-eigenvalues(i)*eye(size(ata,1))).null().mat];
            end
            s = qwe(s);
            v = qwe(v).gram();
            u = sym([]);
            av = this.mat*v.mat;
            for i = 1:min(size(s.mat,1),size(s.mat,2))
                sigma = s.mat(i,i);
                if sigma == 0
                    break
                end
                u = simplify([u av(:,i)/sigma]);
            end
            u = qwe(u);
            u = qwe([u.mat u.lnull().mat]).gram();
        end
        %TODO exact svd

    end
end



