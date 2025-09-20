classdef CompoundOperation
    properties
        oprs
        invs
    end
    methods (Static)
        function [operation] = makeREF(matWrap)
            oprs = {};
            leadingRowNum = 1;
            mat = matWrap.mat;
            for colNum = 1:size(mat,2)
                if leadingRowNum > size(mat,1)
                    break
                end
                entryFound = false;
                for rowNum = leadingRowNum:size(mat,1)
                    val = mat(rowNum,colNum);
                    if val ~= 0
                        if rowNum ~= leadingRowNum
                            mat([rowNum,leadingRowNum],:) = mat([leadingRowNum,rowNum],:);
                            oprs{end+1} = ElementaryOperation.makeRowSwap(rowNum,leadingRowNum);
                        end
                        entryFound = true;
                        break
                    end
                end
                if entryFound & leadingRowNum < size(mat,1)
                    for rowNum = leadingRowNum+1:size(mat,1)
                        if mat(rowNum,colNum) ~= 0
                            scalar = -mat(rowNum,colNum)/mat(leadingRowNum,colNum);
                            mat(rowNum,:) = mat(rowNum,:) + scalar * mat(leadingRowNum,:);
                            oprs{end+1} = ElementaryOperation.makeRowAddition(rowNum,leadingRowNum,scalar);
                        end
                    end
                end
                if entryFound
                    leadingRowNum = leadingRowNum + 1;
                end
            end
            operation = CompoundOperation(oprs);
        end
        function [operation] = makeRREF(matWrap)
            refOperation = CompoundOperation.makeREF(matWrap);
            refMatWrap = matWrap.freeze();
            refOperation.apply(refMatWrap);
            mat = refMatWrap.mat;
            oprs = {};
            %find the rightmost pivot column
            lastPivotColNum = size(mat,2);
            pivotColumnFound = false;   
            for rowNum = size(mat,1):-1:1                
                if pivotColumnFound
                    break
                end
                for colNum = 1:size(mat,2)
                    val = mat(rowNum,colNum);
                    if val ~= 0
                        lastPivotColNum = colNum;
                        pivotColumnFound = true;
                        break
                    end
                end
            end
            for colNum = lastPivotColNum:-1:1
                hasLeadingEntry = false;
                %find row of leading entry and delete preceding rows
                for rowNum = size(mat,1):-1:1
                    val = mat(rowNum,colNum);
                    if hasLeadingEntry
                        scalar = -val;
                        mat(rowNum,:) = mat(rowNum,:) + scalar * mat(leadingEntryRowNum,:);
                        oprs{end+1} = ElementaryOperation.makeRowAddition(rowNum,leadingEntryRowNum,scalar);
                    elseif val ~= 0
                        leadingEntryRowNum = rowNum;
                        scalar = 1/mat(rowNum,colNum);
                        mat(rowNum,:) = mat(rowNum,:) * scalar;
                        oprs{end+1} = ElementaryOperation.makeRowMultiplication(rowNum,scalar);
                        hasLeadingEntry = true;
                    end
                end
            end
            operation = CompoundOperation({refOperation,CompoundOperation(oprs)});
            
        end

    end

    methods
        function obj = CompoundOperation(oprs)
            obj.oprs = oprs;
        end
        function apply(obj, matWrap)
            for i = 1:length(obj.oprs)
                opr = obj.oprs{i};
                opr.apply(matWrap);
            end
        end
        function inverse = invert(obj)
            inverse = {};
            for i = length(obj.oprs):-1:1
                opr = obj.oprs{i};
                inv = opr.invert();
                inverse{end+1} = inv;
            end
        end
        function mop = toMop(obj,n)
            mat = eye(n);
            mop = MatOpsWrapper(mat);
            obj.apply(mop);
        end
    end
end