function [cmap]=redblue()
cmap=[ones(32,1),linspace(0,1,32)',linspace(0,1,32)';linspace(1,0,32)',linspace(1,0,32)',ones(32,1)];