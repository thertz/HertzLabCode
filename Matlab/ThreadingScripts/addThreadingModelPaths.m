function [] = addThreadingModelPaths(threadingModelRoot)

if(~exist('threadingModelRoot','var'))
  threadingModelRoot = 'e:\Work\MyScripts\matlabScripts\ThreadingScripts\';
end

addpath(threadingModelRoot);
addpath([threadingModelRoot,'DataFiles\']);
addpath([threadingModelRoot,'support_functions\']);

