function [ObjectIndexEnv1] = ObjectsIndex(TextQ1Env1,NumQ1Env1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ObjectIndexEnv1 = repmat(0,length(NumQ1Env1),1);
ObjectIndexEnv1(contains(TextQ1Env1,'T1')) = 0;
ObjectIndexEnv1(contains(TextQ1Env1,'F1')) = 1;
ObjectIndexEnv1(contains(TextQ1Env1,'F2')) = 2;
ObjectIndexEnv1(contains(TextQ1Env1,'F3')) = 3;
ObjectIndexEnv1(contains(TextQ1Env1,'F4')) = 4;
ObjectIndexEnv1(contains(TextQ1Env1,'F5')) = 5;
ObjectIndexEnv1(contains(TextQ1Env1,'F6')) = 6;
ObjectIndexEnv1(contains(TextQ1Env1,'F7')) = 7;
ObjectIndexEnv1(contains(TextQ1Env1,'F8')) = 8;
ObjectIndexEnv1(contains(TextQ1Env1,'Catch')) = 9;
end

