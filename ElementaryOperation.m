classdef ElementaryOperation
    properties
        opr
        inv
        oprStr
        invStr
    end
    methods(Static)
        function [operation] = makeRowSwap(r1,r2)
            function opr(matWrap)
                matWrap.mat([r1,r2],:) = matWrap.mat([r2,r1],:);
            end
            oprStr = sprintf("R%d <-> R%d",r1,r2);
            operation = ElementaryOperation(@opr,@opr,oprStr,oprStr);
        end
        function [operation] = makeRowAddition(r1,r2,scalar)
            scalar = simplify(scalar);
            function opr(matWrap)
                matWrap.mat(r1, :) = matWrap.mat(r1, :) + scalar * matWrap.mat(r2, :);
            end
            function inv(matWrap)
                matWrap.mat(r1,:) = matWrap.mat(r1,:) - scalar * matWrap.mat(r2,:);
            end
            oprStr = sprintf("R%d + %s*R%d",r1,string(scalar),r2);
            invStr = sprintf("R%d - %s*R%d",r1,string(scalar),r2);
            operation = ElementaryOperation(@opr, @inv, oprStr, invStr);
        end
        function [operation] = makeRowMultiplication(r1, scalar)
            scalar = simplify(scalar);
            function opr(matWrap)
                matWrap.mat(r1, :) = scalar * matWrap.mat(r1, :);
            end
            function inv(matWrap)
                matWrap.mat(r1, :) = matWrap.mat(r1, :) / scalar;
            end
            oprStr = sprintf("R%d * %s", r1, string(scalar));
            invStr = sprintf("R%d * %s", 1/r1, string(scalar));

         
            operation = ElementaryOperation(@opr, @inv, oprStr, invStr);
        end
    end
    methods
        function obj = ElementaryOperation(opr,inv,oprStr,invStr)
           obj.opr = opr;
           obj.inv = inv;
           obj.oprStr = oprStr;
           obj.invStr = invStr;
        end

        function inverse = invert(obj)
            inverse = ElementaryOperation(obj.inv,obj.opr,obj.invStr,obj.oprStr);
        end

        function apply(obj,matWrap)
            obj.opr(matWrap);
            matWrap.ops{end+1} = obj;
        end
        function str = toString(obj)
            str = obj.oprStr;
        end
        function mop = toMop(obj,n)
            mat = eye(n);
            mop = MatOpsWrapper(mat);
            obj.apply(mop);
        end


    end
end
