%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Copyright 2019 Mohammad Al-Sa'd
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Authors: Mohammad F. Al-Sa'd     (mohammad.al-sad@tuni.fi)
%          Prof. Boualem Boashash  (b.boashash@uq.edu.au)
%
% The following reference should be cited whenever this script is used:
%   M. Al-Sa'd, B. Boashash, "Design and implementation of a multi-sensor
%   newborn EEG seizure and background model with inter-channel field
%   characterization", Digital Signal Processing, 2019.
%
% Last Modification: 25-Jan-2019
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Particle Swarm Optimization
%
% Syntax : out = PSO(params, n, P_real, mP)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% params : Contains the following PSO input parameters:
%          (1) params.MaxIt   : Maximum number of iterations.
%          (2) params.w       : Inertia coefficient.
%          (3) params.wdamp   : Damping ratio of inertia coefficient.
%          (4) params.c1      : Personal acceleration coefficient.
%          (5) params.c2      : Social (Global) acceleration coefficient.
%          (6) params.nPop    : Number of populations.
% n      : number of source signals.
% P_real : Real power map.
% mP     : Mean multi-sensor power.
%
% <OUTPUTs>
% out    : Contains the following PSO outputs:
%          (1) out.BestSol     : Best solution comprised of the following
%                                variables:
%              (1.1) out.BestSol.Position : source signals positions.
%              (1.2) out.BestSol.Cost     : Fitness of the solution.
%          (2) out.GlobalCosts : Global cost at each iteration.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% params.MaxIt = 100;
% params.w     = 0.7298;
% params.wdamp = 0.95;
% params.c1    = 1.49;
% params.c2    = 1.49;
% params.nPop  = 100;
% mP = 200;
% x  = [-2.7 2.7 1.2 ; 0 -3 2 ; 4 0 0];
% P_real = multi_prop_opt(x, mP);
% out = PSO(params, 3, P_real, mP);
% P_sim = multi_prop_opt(out.BestSol.Position, mP);
% figure; plot(P_real); hold on; plot(P_sim,'r');
% legend('real power map','estimated power map');
% figure; plot(out.GlobalCosts);
% xlabel('Iterations'); ylabel('Cost');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function out = PSO(params, n, P_real, mP)
%% Constants
scalp_R = 5.95;              % Head Radius in cm
th1     = 0.29;              % Scalp thickness in cm
th2     = 0.6;               % Skull thickness in cm
th3     = 0.3;               % CSF thickness in cm
brain_R = scalp_R - th1 - th2 - th3;

%% Cost Function
CostFunction = @(x) Cost_Fun(x,P_real, mP);   % Cost function to be minimised

%% Problem Definition
nVar            = 3;                       % Number of unknown (decision) variables
VarSize         = [n nVar];                % Matrix size of decision variables
R_Max           = brain_R-0.01;            % Upper bound of the radius
R_Min           = 0;                       % Lower bound of the radius
MaxVelocity     = 0.2*(R_Max - R_Min);     % Upper bound of velocity
MinVelocity     = -0.2*MaxVelocity;        % Lower bound of velocity

%% Parameters of PSO
MaxIt = params.MaxIt;                % Maximum number of iterations
nPop  = params.nPop;                 % Population size (Swarm size)
w     = params.w;                    % Inertia coefficient
wdamp = params.wdamp;                % Damping ratio of inertia coefficient
c1    = params.c1;                   % Personal acceleration coefficient
c2    = params.c2;                   % Social (Global) acceleration coefficient

%% Initialisation
% Particle template
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost     = [];
% Personal Best
empty_particle.Best.Position = [];
empty_particle.Best.Cost     = [];
% Create Population Array
particle = repmat(empty_particle, nPop, 1);
% Initialise global best
GlobalBest.Cost = Inf;    % Initialised with the worst value. For maximisation problems use -Inf
% Initialise Population Members
for i = 1:nPop
    % Generate random solutions (positions)
    particle(i).Position = event_pos(n);
    % Initialise Velocity
    particle(i).Velocity = zeros(VarSize);
    % Evaluation
    particle(i).Cost = CostFunction(particle(i).Position);
    % Update personal best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    % Update Global Best
    if(particle(i).Best.Cost < GlobalBest.Cost)
        GlobalBest = particle(i).Best;
    end
end
% Array to hold best cost value on each iteration
BestCosts = zeros(MaxIt,1);

%% Main Loop of PSO
for it = 1:MaxIt
    for i = 1:nPop
        % Update particle velocity
        particle(i).Velocity = w*particle(i).Velocity...
            + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position)...
            + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);
        
        % Apply velocity lower and upper bounds limits
        particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
        particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
        % Update particle position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Apply position lower and upper bounds limits
        for nn = 1:n
            [THI,PHI,R] = cart2sph(particle(i).Position(nn,1),particle(i).Position(nn,2),particle(i).Position(nn,3));
            R = min(R,R_Max); R = max(R,R_Min);
            [particle(i).Position(nn,1),particle(i).Position(nn,2),particle(i).Position(nn,3)] = sph2cart(THI,PHI,R);
        end
        % Compute particle cost
        particle(i).Cost = CostFunction(particle(i).Position);
        
        % Update Best Values
        if(particle(i).Cost < particle(i).Best.Cost)
            % Personal best
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            % Global best
            if(particle(i).Best.Cost < GlobalBest.Cost)
                GlobalBest = particle(i).Best;
            end
        end
    end
    
    % Store the best Cost Value
    BestCosts(it) = GlobalBest.Cost;
    % Damping inertia coefficient
    w = w*wdamp;
end
out.BestSol = GlobalBest;
out.GlobalCosts = BestCosts;
end