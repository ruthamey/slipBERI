function [ m_trial ] = do_ensemble_stretch_move( m_curr, walkernumber, otherwalker, z)
% This function takes a current model parameter and performs a 'stretch
% move' to draw a new trial.
% For efficiency, each model parameter has n 'walkers', and each new trial
% is drawn by choosing a value between its current value and another,
% randomly chosen walker. Some stretch beyond the other walker is
% permitted, controlled by parameter a.
% This is based on the method described by Goodman & Weare 2010 in
% 'Communications in Applied Mathematics and Computational Science'
%
% Inputs =
%   m_curr = Current model parameter. Can just be one number, or can be a
%            matrix where each row is another model parameter (e.g. slip on
%            patch 1, slip on patch 2, slip on patch 3...) and each column
%            is another walker
%        a = Value that controls the overreach permitted in a stretch move
%            between one two walkers
% nwalkers = the number of walkers for a model parameter
%
%
% K1 and K2 are parameters for calculating the CDF relating to their
% suggested sampling equation.
% They suggested using PDF for z, p(z)=K/sqrt(z) 1/a>z<a 
% Therefore the CDF, F(z)=K1*(sqrt(z)-1/sqrt(a)), 
% with K1=1/(sqrt(a)-1/sqrt(a))
% so z=(F(z)/K1+K2)^2
%
% Ruth Amey (ruthmjamey@gmail.com) based on a code by Andy Hooper
% 25/8/17
%
% Updated 2/10/17 so z and otherwalker are calcualted at the start of each
% stretch move iteration, to use the same otherwalker for all parameters in
% that stretch move iteration

m_trial(:,1) =m_curr(:,otherwalker)+z.*(m_curr(:,walkernumber)-m_curr(:,otherwalker));

end

