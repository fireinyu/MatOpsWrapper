classdef MatOpsWrapper < handle
    properties
        mat
        ops
        aug = false
    end
    
    methods
        function obj = MatOpsWrapper(mat)
            obj.mat = mat;
            obj.ops = {};
      
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
                separator = NaN(size(obj.mat, 1), 1);
                oprStr = [obj.mat,separator,obj.aug];
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
        end
        function augment(obj,mat)
            obj.aug = MatOpsWrapper(mat);
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
            copyObj = MatOpsWrapper(obj.freezeMat());
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
            mop = oriMop.aug;
        end
        function [l, u] = lu(obj)
            mop = obj.freeze();
            mop.delOps(length(mop.ops));
            mop.ref();
            l = mop.opsMop().invert();
            u = mop;
        end
      
    end
end



