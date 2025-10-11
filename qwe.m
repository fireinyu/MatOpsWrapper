% shortcut class
classdef qwe < MatOpsWrapper
    methods
        function obj = qwe(mat)
            obj@MatOpsWrapper(mat);
        end

        function au(obj,mat)
            obj.augment(mat);
        end

        function res = in(obj)
            res = obj.invert();
        end
        
        function oprStr = rs(obj,r1,r2)
            oprStr = obj.swapRows(r1,r2);
        end
        function oprStr = ra(obj,r1,r2,scalar)
            oprStr = obj.addRow(r1,r2,scalar);
        end
        function oprStr = rm(obj,r1,scalar)
            oprStr = obj.mulRow(r1,scalar);
        end
        function oprStr = re(obj)
            oprStr = obj.ref();
        end
        function oprStr = rr(obj)
            oprStr = obj.rref();
        end

        function sh(obj)
            obj.show();
        end
        function ss(obj)
            obj.steps();
        end
        function si(obj)
            obj.simp();
        end
        function su(obj,value,replacement)
            obj.subs(value,replacement)
        end

        function un(obj,count)
            obj.undo(count);
        end
        function do(obj,upto)
            obj.delOps(upto);
        end

    end
end