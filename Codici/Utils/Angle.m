classdef Angle
% This class is used to define an angle in radians with it's value always
% bound between 0 and 2*pi.
% Definition: Angle(value)
% Properties: 
%       - value: The value of the angle (bounded to [0; 2*pi])
%       - wasNegative: If true, the original value was less than zero
% Methods:
%       - double: Returns the value of the angle as double
%       - tostr:  Returns the value of the angle as string

    properties
        value;
        wasNegative = false;
    end
    methods
        function obj = Angle(Value)
            obj.value = [Value];
            obj = obj.checkRange();
        end
        function obj = checkRange(obj)
            Value = [obj.value];
            for i=1:length(Value)
                if Value(i) < 0
                    obj(i).wasNegative = true;
                    while Value(i)<0
                        Value(i) = Value(i)+2*pi;
                    end
                elseif Value(i) > 2*pi
                    while Value(i) > 2*pi
                        Value(i) = Value(i) - 2*pi;
                    end
                end
            end
            obj.value = Value;
        end
        function num = double(obj)
            num = double([obj.value]);
        end
        function obj = plus(ob1, ob2)
            calc_val = double(ob1) + double(ob2);
            obj = Angle(calc_val);
        end
        function obj = uplus(ob1)
            obj = Angle(double(ob1));
        end
        function obj = minus(ob1, ob2)
            calc_val = double(ob1) - double(ob2);
            obj = Angle(calc_val);
        end
        function obj = uminus(ob1)
            obj = Angle(-double(ob1));
        end
        function num = times(ob1, ob2)
            num = double(ob1) .* double(ob2);
        end
        function obj = mtimes(ob1, ob2)
            calc_val = double(ob1) * double(ob2);
            obj = Angle(calc_val);
        end
        function num = rdivide(ob1, ob2)
            num = double(ob1)/double(ob2);
        end
        function obj = mrdivide(ob1, ob2)
            calc_val = double(ob1)/double(ob2);
            obj = Angle(calc_val);
        end
        function val = lt(ob1, ob2)
            val = double(ob1) < double(ob2);
        end
        function val = gt(ob1, ob2)
            val = double(ob1) > double(ob2);
        end
        function val = le(ob1, ob2)
            val = double(ob1) <= double(ob2);
        end
        function val = ge(ob1, ob2)
            val = double(ob1) >= double(ob2);
        end
        function val = ne(ob1, ob2)
            val = double(ob1) ~= double(ob2);
        end
        function val = eq(ob1, ob2)
            val = double(ob1) == double(ob2);
        end
        function num = sin(ob1)
            num = sin(double(ob1));
        end
        function obj = asin(num)
            obj = asin(double(num));
        end
        function num = cos(ob1)
            num = cos(double(ob1));
        end
        function obj = acos(num)
            obj = acos(double(num));
        end
        function num = tan(ob1)
            num = tan(double(ob1));
        end
        function obj = atan(num)
            obj = atan(double(num));
        end
        function obj = colon(ob1, d, ob2)
            if nargin == 2
                ob2 = Angle(d);
                d = double(1);
            end
            vect = double(ob1):d:double(ob2);
            obj = Angle(vect);
        end
        function obj = linspace(ob1, ob2, num)
            vect = linspace(double(ob1), double(ob2), num);
            obj = Angle(vect);
        end
        function num = size(obj, dim)
            if nargin == 1
                num = length(obj.value);
            else
                num = size(obj.value, dim);
            end
        end
        function obj = transpose(ob1)
            obj = ob1;
            obj.value = ob1.value';
        end
        function obj = ctranspose(ob1)
            obj = ob1;
            obj.value = ob1.value';
        end
        function str = tostr(obj)
            vect = [];
            for i=1:size(obj, 1)
                for j=1:size(obj,2)
                    vect(i,j) = double(obj(i, j));
                end
            end
            str = vect;
            str = string(str);
        end
        function disp(obj)
            disp(obj.value);
        end
    end
end