// Eden Forbes
// MinCogEco Script

// ***************************************
// INCLUDES
// ***************************************

#include "Prey.h"
#include "Predator.h"
#include "random.h"
#include "TSearch.h"
#include <iostream>
#include <iomanip> 
#include <vector>

// ================================================
// A. PARAMETERS & GEN-PHEN MAPPING
// ================================================

// Run constants
// Make sure SpaceSize is also specified in the Prey.cpp and Predator.cpp files
const int SpaceSize = 5000;
const int CoexistTrials = 3;
const int finalCC = 4;
const int CC = 4;
const int minCC = 0;
const int maxCC = 5;
const double G_Rate = 0.01;
// 0-Indexed (0 = 1)
const int start_prey = 20;
const int start_pred = 0;

// Time Constants
// Evolution Step Size:
const double StepSize = 0.1;
// Analysis Step Size:
const double BTStepSize = 0.01;
// Evolution Run Time:
const double RunDuration = 10000;
// Behavioral Trace Run Time:
const double PlotDuration = 10000;
// EcoRate Collection Run Time:
const double RateDuration = 100000;
// Sensory Sample Run Time:
const double SenseDuration = 2000;

// EA params
const int POPSIZE = 100;
const int GENS = 200;
const double MUTVAR = 0.1;
const double CROSSPROB = 0.0;
const double EXPECTED = 1.1;
const double ELITISM = 0.02;

// Nervous system params
const int prey_netsize = 3; 
// weight range 
const double WR = 16.0;
// sensor range
const double SR = 20.0;
// bias range
const double BR = 16.0;
// time constant min max
const double TMIN = 0.5;
const double TMAX = 20.0;
// Weights + TCs & Biases + SensorWeights
const int VectSize = (prey_netsize*prey_netsize + 2*prey_netsize + 3*prey_netsize);

// Prey Sensory Parameters
const double prey_gain = 3.5;
const double prey_s_width = 100.0;

// Prey Metabolic Parameters
// prey_loss_scalar MUST be greater than 1 to prevent strafing behaviors
const double prey_loss_scalar = 3;
const double prey_frate = 0.15;
const double prey_feff = 0.1;
const double prey_repo = 1.5;
const double prey_movecost = 0.0;
const double prey_b_thresh = 3.0;
const double prey_metaloss = ((prey_feff*(finalCC+1))/(SpaceSize/(prey_gain*StepSize))) * prey_loss_scalar;
const double prey_BT_metaloss = ((prey_feff*(finalCC+1))/(SpaceSize/(prey_gain*BTStepSize))) * prey_loss_scalar;

// Predator Sensory Parameters 
const double pred_gain = 3.0;
const double pred_s_width = 110.0;

// Predator Metabolic Parameters
const double pred_frate = 1.0;
const double pred_handling_time = 100.0;
// For Predator Condition:
// 1.0 drift left
// 2.0 drift right
// 3.0 stay still
const double pred_condition = 3.0;

// ------------------------------------
// Genotype-Phenotype Mapping Function
// ------------------------------------

void GenPhenMapping(TVector<double> &gen, TVector<double> &phen)
{
	int k = 1;
	// Prey Time-constants
	for (int i = 1; i <= prey_netsize; i++) {
		phen(k) = MapSearchParameter(gen(k), TMIN, TMAX);
		k++;
	}
	// Prey Bias
	for (int i = 1; i <= prey_netsize; i++) {
		phen(k) = MapSearchParameter(gen(k), -BR, BR);
		k++;
	}
	// Prey Weights
	for (int i = 1; i <= prey_netsize; i++) {
		for (int j = 1; j <= prey_netsize; j++) {
            phen(k) = MapSearchParameter(gen(k), -WR, WR);
			k++;
		}
	}
	// Prey Sensor Weights
	for (int i = 1; i <= 3*prey_netsize; i++) {
        phen(k) = MapSearchParameter(gen(k), -SR, SR);
		k++;
	}
}

// ================================================
// B. TASK ENVIRONMENT & FITNESS FUNCTION
// ================================================

double Coexist(TVector<double> &genotype, RandomState &rs) 
{
    // Set running outcome variable
    double outcome = 99999999999.0;
    // Translate genotype to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);
    // Initialize Prey & Predator agents
    Prey Agent1(prey_netsize, prey_gain, prey_s_width, prey_frate, prey_feff, prey_metaloss, prey_movecost, prey_b_thresh);
    Predator Agent2(pred_gain, pred_s_width, pred_frate, pred_handling_time);
    // Set Prey nervous system
    Agent1.NervousSystem.SetCircuitSize(prey_netsize);
    int k = 1;
    // Prey Time-constants
    for (int i = 1; i <= prey_netsize; i++) {
        Agent1.NervousSystem.SetNeuronTimeConstant(i,phenotype(k));
        k++;
    }
    // Prey Biases
    for (int i = 1; i <= prey_netsize; i++) {
        Agent1.NervousSystem.SetNeuronBias(i,phenotype(k));
        k++;
    }
    // Prey Neural Weights
    for (int i = 1; i <= prey_netsize; i++) {
        for (int j = 1; j <= prey_netsize; j++) {
            Agent1.NervousSystem.SetConnectionWeight(i,j,phenotype(k));
            k++;
        }
    }
    int j = k;
    // Prey Sensor Weights
    for (int i = 1; i <= prey_netsize*3; i++) {
        Agent1.sensorweights[i] = phenotype(k);
        k++;
    }
    // Run Simulation
    for (int trial = 0; trial < 3*CoexistTrials; trial++){
        // Set Predator Condition for Trial
        double pred_condition = 0.0;
        if(trial < CoexistTrials){
            pred_condition = 1.0;
        }
        else if(trial >= CoexistTrials && trial < 2*CoexistTrials){
            pred_condition = 2.0;
        }
        else if(trial >= 2*CoexistTrials && trial < 3*CoexistTrials){
            pred_condition = 3.0;
        }
        Agent2.condition = pred_condition;
        // Reset Prey agent, randomize its location
        Agent1.Reset(rs.UniformRandomInteger(0,SpaceSize), 1.0);
        // Seed preylist with starting population
        TVector<Prey> preylist(0,0);
        preylist[0] = Agent1;
        for (int i = 0; i < start_prey; i++){
            Prey newprey(prey_netsize, prey_gain, prey_s_width, prey_frate, prey_feff, prey_metaloss, prey_movecost, prey_b_thresh);
            // Reset Prey agent, randomize its location
            newprey.Reset(rs.UniformRandomInteger(0,SpaceSize), 1.0);
            // Copy over nervous system
            newprey.NervousSystem = Agent1.NervousSystem;
            newprey.sensorweights = Agent1.sensorweights;
            // Add to preylist
            preylist.SetBounds(0, preylist.Size());
            preylist[preylist.Size()-1] = newprey;
        }
        // Seed predlist with starting population
        TVector<Predator> predlist(0,0);
        predlist[0] = Agent2;
        for (int i = 0; i < start_pred; i++){
            Predator newpred(pred_gain, pred_s_width, pred_frate, pred_handling_time);
            // Reset Predator agent, randomize its location
            newpred.Reset(rs.UniformRandomInteger(0,SpaceSize));
            newpred.condition = pred_condition;
            // Add to predlist
            predlist.SetBounds(0, predlist.Size());
            predlist[predlist.Size()-1] = newpred;
        }
        // Initialize Producers, fill world to carrying capacity
        TVector<double> food_pos;
        TVector<double> WorldFood(1, SpaceSize);
        WorldFood.FillContents(0.0);
        for (int i = 0; i <= CC; i++){
            int f = rs.UniformRandomInteger(1,SpaceSize);
            WorldFood[f] = 1.0;
            food_pos.SetBounds(0, food_pos.Size());
            food_pos[food_pos.Size()-1] = f;
        }
        // Set Clocks & trial outcome variables
        double clock = 0.0;
        double prey_snacks = 0.0;
        double pred_snacks = 0.0;
        double prey_outcome = RunDuration;
        double pred_outcome = RunDuration;
        double running_pop = 0.0;

        // Run a Trial
        for (double time = 0; time < RunDuration; time += StepSize){
            // Remove any consumed food from food list
            TVector<double> dead_food(0,-1);
            for (int i = 0; i < food_pos.Size(); i++){
                if (WorldFood[food_pos[i]] <= 0){
                    dead_food.SetBounds(0, dead_food.Size());
                    dead_food[dead_food.Size()-1] = food_pos[i];
                }
            }
            if (dead_food.Size() > 0){
                for (int i = 0; i < dead_food.Size(); i++){
                    food_pos.RemoveFood(dead_food[i]);
                    food_pos.SetBounds(0, food_pos.Size()-2);
                }
            }
            // Chance for new food to grow
            // Carrying capacity is 0 indexed, add 1 for true amount
            double food_count = food_pos.Size();
            double s_chance = food_count/(CC+1);
            double c = rs.UniformRandom(0,1);
            if (c > s_chance){
                int f = rs.UniformRandomInteger(1,SpaceSize);
                WorldFood[f] = 1.0;
                food_pos.SetBounds(0, food_pos.Size());
                food_pos[food_pos.Size()-1] = f;
            }
            // Update Prey Positions
            TVector<double> prey_pos;
            for (int i = 0; i < preylist.Size(); i++){
                prey_pos.SetBounds(0, prey_pos.Size());
                prey_pos[prey_pos.Size()-1] = preylist[i].pos;
            }
            // Predator Sense & Step
            TVector<Predator> newpredlist;
            TVector<int> preddeaths;
            for (int i = 0; i < predlist.Size(); i++){
                predlist[i].Sense(prey_pos);
                predlist[i].Step(StepSize, WorldFood, preylist);
                if (predlist[i].snackflag > 0){
                    pred_snacks += predlist[i].snackflag;
                    predlist[i].snackflag = 0;
                }
            }
            // Update Predator Positions
            TVector<double> pred_pos;
            for (int i = 0; i < predlist.Size(); i++){
                pred_pos.SetBounds(0, pred_pos.Size());
                pred_pos[pred_pos.Size()-1] = predlist[i].pos;
            }
            // Prey Sense & Step
            TVector<Prey> newpreylist;
            TVector<int> preydeaths;
            for (int i = 0; i < preylist.Size(); i++){
                preylist[i].Sense(food_pos, pred_pos);
                preylist[i].Step(StepSize, WorldFood);
                if (preylist[i].snackflag > 0){
                    prey_snacks += preylist[i].snackflag;
                    preylist[i].snackflag = 0;
                }
                if (preylist[i].birth == true){
                    preylist[i].state = preylist[i].state - prey_repo;
                    preylist[i].birth = false;
                    Prey newprey(prey_netsize, prey_gain, prey_s_width, prey_frate, prey_feff, prey_metaloss, prey_movecost, prey_b_thresh);
                    newprey.NervousSystem = preylist[i].NervousSystem;
                    newprey.sensorweights = preylist[i].sensorweights;
                    newprey.Reset(preylist[i].pos+2, prey_repo);
                    newpreylist.SetBounds(0, newpreylist.Size());
                    newpreylist[newpreylist.Size()-1] = newprey;
                }
                if (preylist[i].death == true){
                    preydeaths.SetBounds(0, preydeaths.Size());
                    preydeaths[preydeaths.Size()-1] = i;
                }
            }
            // Update clocks
            clock += StepSize;
            // Update prey list with new prey list and deaths
            if (preydeaths.Size() > 0){
                for (int i = 0; i <= preydeaths.Size()-1; i++){
                    preylist.RemoveItem(preydeaths[i]);
                    preylist.SetBounds(0, preylist.Size()-2);
                }
            }
            if (newpreylist.Size() > 0){
                for (int i = 0; i <= newpreylist.Size()-1; i++){
                    preylist.SetBounds(0, preylist.Size());
                    preylist[preylist.Size()-1] = newpreylist[i];
                }
            }
            // Check for prey population collapse
            if (preylist.Size() <= 0){
                break;
            }
            // Reset lists for next step
            else{
                running_pop += preylist.Size();
                newpreylist.~TVector();
                preydeaths.~TVector();
                newpredlist.~TVector();
                preddeaths.~TVector();
                prey_pos.~TVector();
                pred_pos.~TVector();
                dead_food.~TVector();
            }
        }
        // Fitness part 1 is proportion of RunTime survived
        double runmeasure = (clock/RunDuration);
        // Fitness part 2 is average population across the run
        double popmeasure = (running_pop/(RunDuration/StepSize));
        // Generate fitness
        double fitmeasure = runmeasure + popmeasure/10;
        // Keep minimum fitness value across trials
        if (fitmeasure < outcome){
            outcome = fitmeasure;
        }
    }
    return outcome;
}

// ================================================
// C. ADDITIONAL EVOLUTIONARY FUNCTIONS
// ================================================
int CCTerminationFunction(int Generation, double BestPerf, double AvgPerf, double PerfVar)
{
	if (BestPerf >= 1.0) return 1;
	else return 0;
}

int EndTerminationFunction(int Generation, double BestPerf, double AvgPerf, double PerfVar)
{
	if (BestPerf >= 100.0) return 1;
	else return 0;
}

void EvolutionaryRunDisplay(int Generation, double BestPerf, double AvgPerf, double PerfVar)
{
	cout << Generation << " " << BestPerf << " " << AvgPerf << " " << PerfVar << endl;
}

void ResultsDisplay(TSearch &s)
{
	TVector<double> bestVector;
	ofstream BestIndividualFile;

	// Save the genotype of the best individual
	bestVector = s.BestIndividual();
	BestIndividualFile.open("best.gen.dat");
    BestIndividualFile << setprecision(32);
	BestIndividualFile << bestVector << endl;
	BestIndividualFile.close();
}

// ================================================
// D. ANALYSIS FUNCTIONS
// ================================================
// ------------------------------------
// Behavioral Traces
// ------------------------------------

double BehavioralTracesCoexist (TVector<double> &genotype, RandomState &rs, double condition) 
{
    // Start output files
	ofstream preyfile("analysis_results/prey_pos.dat");
    ofstream predfile("analysis_results/pred_pos.dat");
    ofstream preypopfile("analysis_results/prey_pop.dat");
    ofstream predpopfile("analysis_results/pred_pop.dat");
	ofstream foodfile("analysis_results/food_pos.dat");
    ofstream foodpopfile("analysis_results/food_pop.dat");
    ofstream bestphen("best.phen.dat");
    // Set running outcome
    double outcome = 99999999999.0;
    // Translate to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);
    // Create agents
    Prey Agent1(prey_netsize, prey_gain, prey_s_width, prey_frate, prey_feff, prey_BT_metaloss, prey_movecost, prey_b_thresh);
    Predator Agent2(pred_gain, pred_s_width, pred_frate,pred_handling_time);
    // Set nervous system
    Agent1.NervousSystem.SetCircuitSize(prey_netsize);
    int k = 1;
    // Prey Time-constants
    for (int i = 1; i <= prey_netsize; i++) {
        Agent1.NervousSystem.SetNeuronTimeConstant(i,phenotype(k));
        k++;
    }
    // Prey Biases
    for (int i = 1; i <= prey_netsize; i++) {
        Agent1.NervousSystem.SetNeuronBias(i,phenotype(k));
        k++;
    }
    // Prey Neural Weights
    for (int i = 1; i <= prey_netsize; i++) {
        for (int j = 1; j <= prey_netsize; j++) {
            Agent1.NervousSystem.SetConnectionWeight(i,j,phenotype(k));
            k++;
        }
    }
    // Prey Sensor Weights
    for (int i = 1; i <= prey_netsize*3; i++) {
        Agent1.sensorweights[i] = phenotype(k);
        k++;
    }
    // Run Simulation
    double prey_outcome = 99999999999.0;
    double pred_outcome = 99999999999.0;
    Agent2.condition = condition;
    // Reset Agents & Vectors
    Agent1.Reset(rs.UniformRandomInteger(0,SpaceSize), 1.5);
    Agent2.Reset(rs.UniformRandomInteger(0,SpaceSize));
    // Seed preylist with starting population
    TVector<Prey> preylist(0,0);
    preylist[0] = Agent1;
    for (int i = 0; i < start_prey; i++){
        Prey newprey(prey_netsize, prey_gain, prey_s_width, prey_frate, prey_feff, prey_metaloss, prey_movecost, prey_b_thresh);
        newprey.Reset(rs.UniformRandomInteger(0,SpaceSize), 1.0);
        newprey.NervousSystem = Agent1.NervousSystem;
        newprey.sensorweights = Agent1.sensorweights;
        preylist.SetBounds(0, preylist.Size());
        preylist[preylist.Size()-1] = newprey;
    }
    // Seed predlist with starting population
    TVector<Predator> predlist(0,0);
    predlist[0] = Agent2;
    for (int i = 0; i < start_pred; i++){
        Predator newpred(pred_gain, pred_s_width, pred_frate, pred_handling_time);
        newpred.Reset(rs.UniformRandomInteger(0,SpaceSize));
        newpred.condition = pred_condition;
        predlist.SetBounds(0, predlist.Size());
        predlist[predlist.Size()-1] = newpred;
    }
    // Fill World to Carrying Capacity
    TVector<double> food_pos(0,-1);
    TVector<double> WorldFood(1, SpaceSize);
    WorldFood.FillContents(0.0);
    for (int i = 0; i <= CC; i++){
        int f = rs.UniformRandomInteger(1,SpaceSize);
        WorldFood[f] = 1.0;
        food_pos.SetBounds(0, food_pos.Size());
        food_pos[food_pos.Size()-1] = f;
    }
    // Run Simulation
    for (double time = 0; time < PlotDuration; time += BTStepSize){
        // Remove chomped food from food list
        TVector<double> dead_food(0,-1);
        for (int i = 0; i < food_pos.Size(); i++){
            if (WorldFood[food_pos[i]] <= 0){
                dead_food.SetBounds(0, dead_food.Size());
                dead_food[dead_food.Size()-1] = food_pos[i];
            }
        }
        if (dead_food.Size() > 0){
            for (int i = 0; i < dead_food.Size(); i++){
                food_pos.RemoveFood(dead_food[i]);
                food_pos.SetBounds(0, food_pos.Size()-2);
            }
        }
        // Carrying capacity is 0 indexed, add 1 for true amount
        double food_count = food_pos.Size();
        double s_chance = food_count/(CC+1);
        double c = rs.UniformRandom(0,1);
        if (c > s_chance){
            int f = rs.UniformRandomInteger(1,SpaceSize);
            WorldFood[f] = 1.0;
            food_pos.SetBounds(0, food_pos.Size());
            food_pos[food_pos.Size()-1] = f;
        }
        // Update Prey Positions
        TVector<double> prey_pos;
        for (int i = 0; i < preylist.Size(); i++){
            prey_pos.SetBounds(0, prey_pos.Size());
            prey_pos[prey_pos.Size()-1] = preylist[i].pos;
        }
        // Predator Sense & Step
        TVector<Predator> newpredlist;
        TVector<int> preddeaths;
        for (int i = 0; i < predlist.Size(); i++){
            predlist[i].Sense(prey_pos);
            predlist[i].Step(BTStepSize, WorldFood, preylist);
        }
        // Update Predator Positions
        TVector<double> pred_pos;
        for (int i = 0; i < predlist.Size(); i++){
            pred_pos.SetBounds(0, pred_pos.Size());
            pred_pos[pred_pos.Size()-1] = predlist[i].pos;
        }
        // Prey Sense & Step
        TVector<Prey> newpreylist;
        TVector<int> preydeaths;
        for (int i = 0; i < preylist.Size(); i++){
            preylist[i].Sense(food_pos, pred_pos);
            preylist[i].Step(BTStepSize, WorldFood);
            if (preylist[i].birth == true){
                preylist[i].state = preylist[i].state - prey_repo;
                preylist[i].birth = false;
                Prey newprey(prey_netsize, prey_gain, prey_s_width, prey_frate, prey_feff, prey_metaloss, prey_movecost, prey_b_thresh);
                newprey.NervousSystem = preylist[i].NervousSystem;
                newprey.sensorweights = preylist[i].sensorweights;
                newprey.Reset(preylist[i].pos+2, prey_repo);
                newpreylist.SetBounds(0, newpreylist.Size());
                newpreylist[newpreylist.Size()-1] = newprey;
            }
            if (preylist[i].death == true){
                preydeaths.SetBounds(0, preydeaths.Size());
                preydeaths[preydeaths.Size()-1] = i;
            }
        }
        // Update prey list with new prey list and deaths
        if (preydeaths.Size() > 0){
            for (int i = 0; i < preydeaths.Size(); i++){
                preylist.RemoveItem(preydeaths[i]);
                preylist.SetBounds(0, preylist.Size()-2);
            }
        }
        if (newpreylist.Size() > 0){
            for (int i = 0; i < newpreylist.Size(); i++){
                preylist.SetBounds(0, preylist.Size());
                preylist[preylist.Size()-1] = newpreylist[i];
            }
        }
        // Save
        preyfile << prey_pos << endl;
        preypopfile << preylist.Size() << " ";
        predfile << pred_pos << endl;
        predpopfile << predlist.Size() << " ";
        foodfile << food_pos << endl;
        double foodsum = 0.0;
        for (int i = 0; i < food_pos.Size(); i++){
            foodsum += WorldFood[food_pos[i]];
        }
        foodpopfile << foodsum << " ";
        // Check Population Collapse
        if (preylist.Size() <= 0){
            break;
        }
        else{
            newpreylist.~TVector();
            preydeaths.~TVector();
            newpredlist.~TVector();
            preddeaths.~TVector();
            prey_pos.~TVector();
            pred_pos.~TVector();
            dead_food.~TVector();
        }
    }
    preyfile.close();
    preypopfile.close();
    predfile.close();
    predpopfile.close();
	foodfile.close();
    foodpopfile.close();
    // Save best phenotype
    bestphen << phenotype << endl;
    return 0;
}

// ------------------------------------
// Interaction Rate Data Functions
// ------------------------------------
void DeriveLambdaH(Prey &prey, Predator &predator, RandomState &rs, double &maxCC, int &samplesize, double &transient)
{
    ofstream lambHfile("analysis_results/lambH.dat");
    for (int j = -1; j <= maxCC; j++){
        TVector<double> lambH;
        for (int k = 0; k <= samplesize; k++){
            int carrycapacity = j;
            // Fill World to Carrying Capacity
            TVector<double> food_pos;
            TVector<double> WorldFood(1, SpaceSize);
            WorldFood.FillContents(0.0);
            for (int i = 0; i <= CC; i++){
                int f = rs.UniformRandomInteger(1,SpaceSize);
                WorldFood[f] = 1.0;
                food_pos.SetBounds(0, food_pos.Size());
                food_pos[food_pos.Size()-1] = f;
            }
            // Make dummy predator list
            TVector<double> pred_pos(0,-1);
            double munch_count = 0;
            for (double time = 0; time < RateDuration; time += StepSize){
                // Remove chomped food from food list
                TVector<double> dead_food(0,-1);
                for (int i = 0; i < food_pos.Size(); i++){
                    if (WorldFood[food_pos[i]] <= 0){
                        dead_food.SetBounds(0, dead_food.Size());
                        dead_food[dead_food.Size()-1] = food_pos[i];
                    }
                }
                if (dead_food.Size() > 0){
                    for (int i = 0; i < dead_food.Size(); i++){
                        food_pos.RemoveFood(dead_food[i]);
                        food_pos.SetBounds(0, food_pos.Size()-2);
                    }
                }
                // Carrying capacity is 0 indexed, add 1 for true amount
                double food_count = food_pos.Size();
                double s_chance = G_Rate*food_count/(CC+1);
                double c = rs.UniformRandom(0,1);
                if (c > s_chance){
                    int f = rs.UniformRandomInteger(1,SpaceSize);
                    WorldFood[f] = 1.0;
                    food_pos.SetBounds(0, food_pos.Size());
                    food_pos[food_pos.Size()-1] = f;
                }
                // Prey Sense & Step
                prey.Sense(food_pos, pred_pos);
                prey.Step(StepSize, WorldFood);
                // Check Births
                if (prey.birth == true){
                    prey.state = prey.state - prey_repo;
                    prey.birth = false;
                }
                // Check Deaths
                if (prey.death == true){
                    prey.Reset(rs.UniformRandomInteger(0,SpaceSize), 1.5);
                    prey.death = false;
                }
                // Check # of times food crossed
                if (time > transient){
                    munch_count += prey.snackflag;
                    prey.snackflag = 0.0;
                }
            }
            double munchrate = munch_count/((RateDuration-transient)/StepSize);
            lambH.SetBounds(0, lambH.Size());
            lambH[lambH.Size()-1] = munchrate;
        }
        lambHfile << lambH << endl;
        lambH.~TVector();
    }
    // Save
    lambHfile.close();
}

void DeriveLambdaP(Prey &prey, Predator &predator, RandomState &rs, double &maxprey, int &samplesize, double &transient)
{
    ofstream lambCfile("analysis_results/lambC.dat");
    for (int j = -1; j<=maxprey; j++)
    {   
        TVector<double> lambC;
        for (int k = 0; k<=samplesize; k++){
            // Fill World to Carrying Capacity
            TVector<double> food_pos;
            TVector<double> WorldFood(1, SpaceSize);
            WorldFood.FillContents(0.0);
            for (int i = 0; i <= CC; i++){
                int f = rs.UniformRandomInteger(1,SpaceSize);
                WorldFood[f] = 1.0;
                food_pos.SetBounds(0, food_pos.Size());
                food_pos[food_pos.Size()-1] = f;
            }
            // Seed preylist with starting population
            TVector<Prey> preylist(0,0);
            TVector<double> prey_pos;
            preylist[0] = prey;
            for (int i = 0; i < j; i++){
                Prey newprey(prey_netsize, prey_gain, prey_s_width, prey_frate, prey_feff, prey_metaloss, prey_movecost, prey_b_thresh);
                newprey.Reset(rs.UniformRandomInteger(0,SpaceSize), 1.5);
                newprey.NervousSystem = prey.NervousSystem;
                newprey.sensorweights = prey.sensorweights;
                preylist.SetBounds(0, preylist.Size());
                preylist[preylist.Size()-1] = newprey;
                }
            double munch_count = 0;
            for (double time = 0; time < RateDuration; time += StepSize){
                // Remove chomped food from food list
                TVector<double> dead_food(0,-1);
                for (int i = 0; i < food_pos.Size(); i++){
                    if (WorldFood[food_pos[i]] <= 0){
                        dead_food.SetBounds(0, dead_food.Size());
                        dead_food[dead_food.Size()-1] = food_pos[i];
                    }
                }
                if (dead_food.Size() > 0){
                    for (int i = 0; i < dead_food.Size(); i++){
                        food_pos.RemoveFood(dead_food[i]);
                        food_pos.SetBounds(0, food_pos.Size()-2);
                    }
                }
                // Carrying capacity is 0 indexed, add 1 for true amount
                double food_count = food_pos.Size();
                double s_chance = G_Rate*food_count/(CC+1);
                double c = rs.UniformRandom(0,1);
                if (c > s_chance){
                    int f = rs.UniformRandomInteger(1,SpaceSize);
                    WorldFood[f] = 1.0;
                    food_pos.SetBounds(0, food_pos.Size());
                    food_pos[food_pos.Size()-1] = f;
                }
                // Prey Sense & Step
                TVector<Prey> newpreylist;
                TVector<int> deaths;
                TVector<double> prey_pos;
                TVector<double> pred_pos;
                pred_pos.SetBounds(0, pred_pos.Size());
                pred_pos[pred_pos.Size()-1] = predator.pos;
                for (int i = 0; i < preylist.Size(); i++){
                    if (preylist[i].death == true){
                        preylist[i].Reset(rs.UniformRandomInteger(0,SpaceSize), 1.5);
                    }
                    else{
                        preylist[i].Sense(food_pos, pred_pos);
                        preylist[i].Step(StepSize, WorldFood);
                    }
                }
                for (int i = 0; i <= preylist.Size()-1; i++){
                    prey_pos.SetBounds(0, prey_pos.Size());
                    prey_pos[prey_pos.Size()-1] = preylist[i].pos;
                }

                // Predator Sense & Step
                predator.Sense(prey_pos);
                predator.Step(StepSize, WorldFood, preylist);
                // Check # of times food crossed
                if(time > transient){
                    munch_count += predator.snackflag;
                    predator.snackflag = 0.0;
                }
            }

            double munchrate = munch_count/((RateDuration-transient)/StepSize);
            lambC.SetBounds(0, lambC.Size());
            lambC[lambC.Size()-1] = munchrate;
        }
        lambCfile << lambC << endl;
        lambC.~TVector();
    }
    // Save
    lambCfile.close();
}

void CollectEcoRates(TVector<double> &genotype, RandomState &rs)
{
    ofstream erates("analysis_results/ecosystem_rates.dat");
    // Translate to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);
    // Create agents
    Prey Agent1(prey_netsize, prey_gain, prey_s_width, prey_frate, prey_feff, prey_metaloss, prey_movecost, prey_b_thresh);
    Predator Agent2(pred_gain, pred_s_width, pred_frate, pred_handling_time);
    Agent2.condition = pred_condition;
    // Set nervous system
    Agent1.NervousSystem.SetCircuitSize(prey_netsize);
    int k = 1;
    // Prey Time-constants
    for (int i = 1; i <= prey_netsize; i++) {
        Agent1.NervousSystem.SetNeuronTimeConstant(i,phenotype(k));
        k++;
    }
    // Prey Biases
    for (int i = 1; i <= prey_netsize; i++) {
        Agent1.NervousSystem.SetNeuronBias(i,phenotype(k));
        k++;
    }
    // Prey Neural Weights
    for (int i = 1; i <= prey_netsize; i++) {
        for (int j = 1; j <= prey_netsize; j++) {
            Agent1.NervousSystem.SetConnectionWeight(i,j,phenotype(k));
            k++;
        }
    }
    // Prey Sensor Weights
    for (int i = 1; i <= prey_netsize*3; i++) {
        Agent1.sensorweights[i] = phenotype(k);
        k++;
    }
    // Save Growth Rates
    // Max growth rate of producers is the chance of a new plant coming in on a given time step
    double systemcc = CC+1; // 0 indexed
    double rr = 1;
    erates << rr << " ";
    erates << systemcc << " ";
    erates << Agent1.frate << " ";
    erates << Agent1.feff << " ";
    erates << Agent1.metaloss << " ";
    erates << Agent2.frate << " ";
    // erates << Agent2.feff << " ";
    // erates << Agent2.metaloss << " ";
    erates.close();

    // Set Sampling Range & Frequency
    double maxCC = 50;
    double maxprey = 50;
    double transient = 5000.0;
    int samplesize = 6;

    // Collect Prey Lambda & r
    printf("Collecting Prey rates\n");
    DeriveLambdaH(Agent1, Agent2, rs, maxCC, samplesize, transient);
    // Collect Predator Lambda & r
    printf("Collecting Predator rates\n");
    DeriveLambdaP(Agent1, Agent2, rs, maxprey, samplesize, transient);
}

// ------------------------------------
// Sensory Sample Functions
// ------------------------------------
void SensorySample(TVector<double> &genotype, RandomState &rs)
{
    // Translate to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);
    // Create agents
    Prey Agent1(prey_netsize, prey_gain, prey_s_width, prey_frate, prey_feff, prey_BT_metaloss, prey_movecost, prey_b_thresh);
    Predator Agent2(pred_gain, pred_s_width, pred_frate, pred_handling_time);
    Agent2.condition = pred_condition;
    // Set nervous system
    Agent1.NervousSystem.SetCircuitSize(prey_netsize);
    int k = 1;
    // Prey Time-constants
    for (int i = 1; i <= prey_netsize; i++) {
        Agent1.NervousSystem.SetNeuronTimeConstant(i,phenotype(k));
        k++;
    }
    // Prey Biases
    for (int i = 1; i <= prey_netsize; i++) {
        Agent1.NervousSystem.SetNeuronBias(i,phenotype(k));
        k++;
    }
    // Prey Neural Weights
    for (int i = 1; i <= prey_netsize; i++) {
        for (int j = 1; j <= prey_netsize; j++) {
            Agent1.NervousSystem.SetConnectionWeight(i,j,phenotype(k));
            k++;
        }
    }
    // Prey Sensor Weights
    for (int i = 1; i <= prey_netsize*3; i++) {
        Agent1.sensorweights[i] = phenotype(k);
        k++;
    }
    Agent1.Reset(2000, 1.5);
    Agent2.condition = pred_condition;
    ofstream SSfile("analysis_results/SenS.dat");
    // Fill World to Carrying Capacity
    TVector<double> food_pos(0,-1);
    TVector<double> WorldFood(1, SpaceSize);
    WorldFood.FillContents(0.0);
    for (int i = 0; i <= CC; i++){
        int f = rs.UniformRandomInteger(1,SpaceSize);
        WorldFood[f] = 1.0;
        food_pos.SetBounds(0, food_pos.Size());
        food_pos[food_pos.Size()-1] = f;
    }
    double munch_count = 0;
    TVector<Prey> preylist(0,0);
    preylist[0] = Agent1;
    TVector<double> FS;
    TVector<double> PS;
    TVector<double> SS;
    TVector<double> NO1;
    TVector<double> N1FS;
    TVector<double> N1PS;
    TVector<double> N1SS;
    TVector<double> NO2;
    TVector<double> N2FS;
    TVector<double> N2PS;
    TVector<double> N2SS;
    TVector<double> NO3;
    TVector<double> N3FS;
    TVector<double> N3PS;
    TVector<double> N3SS;
    TVector<double> mov;

    // Simulation 1: Food Sense Sample
    for (double time = 0; time < SenseDuration; time += BTStepSize){
        TVector<double> dead_food(0,-1);
        for (int i = 0; i < food_pos.Size(); i++){
            if (WorldFood[food_pos[i]] <= 0){
                dead_food.SetBounds(0, dead_food.Size());
                dead_food[dead_food.Size()-1] = food_pos[i];
            }
        }
        if (dead_food.Size() > 0){
            for (int i = 0; i < dead_food.Size(); i++){
                food_pos.RemoveFood(dead_food[i]);
                food_pos.SetBounds(0, food_pos.Size()-2);
            }
        }
        double food_count = food_pos.Size();
        double s_chance = food_count/(1);
        double c = rs.UniformRandom(0,1);
        if (c > s_chance){
            int f = rs.UniformRandomInteger(1,SpaceSize);
            WorldFood[f] = 1.0;
            food_pos.SetBounds(0, food_pos.Size());
            food_pos[food_pos.Size()-1] = f;
        }
        // Carrying capacity is 0 indexed, add 1 for true amount
        TVector<double> pred_pos(0,-1);
        Agent1.Sense(food_pos, pred_pos);
        FS.SetBounds(0, SS.Size());
        FS[FS.Size()-1] = Agent1.f_sensor;
        N1FS.SetBounds(0, N1FS.Size());
        N1FS[N1FS.Size()-1] = Agent1.f_sensor * Agent1.sensorweights[1];
        N2FS.SetBounds(0, N2FS.Size());
        N2FS[N2FS.Size()-1] = Agent1.f_sensor * Agent1.sensorweights[4];
        N3FS.SetBounds(0, N3FS.Size());
        N3FS[N3FS.Size()-1] = Agent1.f_sensor * Agent1.sensorweights[7];
        PS.SetBounds(0, PS.Size());
        PS[PS.Size()-1] = Agent1.p_sensor;
        N1PS.SetBounds(0, N1PS.Size());
        N1PS[N1PS.Size()-1] = Agent1.p_sensor * Agent1.sensorweights[2];
        N2PS.SetBounds(0, N2PS.Size());
        N2PS[N2PS.Size()-1] = Agent1.p_sensor * Agent1.sensorweights[5];
        N3PS.SetBounds(0, N3PS.Size());
        N3PS[N3PS.Size()-1] = Agent1.p_sensor * Agent1.sensorweights[8];
        SS.SetBounds(0, SS.Size());
        SS[SS.Size()-1] = Agent1.state;
        N1SS.SetBounds(0, N1SS.Size());
        N1SS[N1SS.Size()-1] = Agent1.state * Agent1.sensorweights[3];
        N2SS.SetBounds(0, N2SS.Size());
        N2SS[N2SS.Size()-1] = Agent1.state * Agent1.sensorweights[6];
        N3SS.SetBounds(0, N3SS.Size());
        N3SS[N3SS.Size()-1] = Agent1.state * Agent1.sensorweights[9];

        Agent1.Step(BTStepSize, WorldFood);
        NO1.SetBounds(0, NO1.Size());
        NO1[NO1.Size()-1] = Agent1.NervousSystem.NeuronOutput(1);
        NO2.SetBounds(0, NO2.Size());
        NO2[NO2.Size()-1] = Agent1.NervousSystem.NeuronOutput(2);
        NO3.SetBounds(0, NO3.Size());
        NO3[NO3.Size()-1] = Agent1.NervousSystem.NeuronOutput(3);
        mov.SetBounds(0, mov.Size());
        mov[mov.Size()-1] = (Agent1.NervousSystem.NeuronOutput(2) - Agent1.NervousSystem.NeuronOutput(1));

        if (Agent1.birth == true){
            Agent1.state = Agent1.state - prey_repo;
            Agent1.birth = false;
        }
        if (Agent1.death == true){
            Agent1.Reset(rs.UniformRandomInteger(0,SpaceSize), 1.5);
            Agent1.death = false;
        }

        pred_pos.~TVector();
        dead_food.~TVector();
    }
    WorldFood.FillContents(0.0);
    Agent1.Reset(100, 1.5);
    Agent2.Reset(250);

    // Simulation 2: Predator Sense Sample
    for (double time = 0; time < SenseDuration; time += BTStepSize){
        TVector<double> food_pos(0,-1);
        TVector<double> dead_food(0,-1);
        for (int i = 0; i < food_pos.Size(); i++){
            if (WorldFood[food_pos[i]] <= 0){
                dead_food.SetBounds(0, dead_food.Size());
                dead_food[dead_food.Size()-1] = food_pos[i];
            }
        }
        if (dead_food.Size() > 0){
            for (int i = 0; i < dead_food.Size(); i++){
                food_pos.RemoveFood(dead_food[i]);
                food_pos.SetBounds(0, food_pos.Size()-2);
            }
        }
        // Update Prey Positions
        TVector<double> prey_pos;
        for (int i = 0; i < preylist.Size(); i++){
            prey_pos.SetBounds(0, prey_pos.Size());
            prey_pos[prey_pos.Size()-1] = preylist[i].pos;
        }
        Agent2.Sense(prey_pos);
        Agent2.Step(BTStepSize, WorldFood, preylist);

        TVector<double> pred_pos(0,-1);
        pred_pos.SetBounds(0, pred_pos.Size());
        pred_pos[pred_pos.Size()-1] = Agent2.pos;

        Agent1.Sense(food_pos, pred_pos);
        FS.SetBounds(0, SS.Size());
        FS[FS.Size()-1] = Agent1.f_sensor;
        N1FS.SetBounds(0, N1FS.Size());
        N1FS[N1FS.Size()-1] = Agent1.f_sensor * Agent1.sensorweights[1];
        N2FS.SetBounds(0, N2FS.Size());
        N2FS[N2FS.Size()-1] = Agent1.f_sensor * Agent1.sensorweights[4];
        N3FS.SetBounds(0, N3FS.Size());
        N3FS[N3FS.Size()-1] = Agent1.f_sensor * Agent1.sensorweights[7];
        PS.SetBounds(0, PS.Size());
        PS[PS.Size()-1] = Agent1.p_sensor;
        N1PS.SetBounds(0, N1PS.Size());
        N1PS[N1PS.Size()-1] = Agent1.p_sensor * Agent1.sensorweights[2];
        N2PS.SetBounds(0, N2PS.Size());
        N2PS[N2PS.Size()-1] = Agent1.p_sensor * Agent1.sensorweights[5];
        N3PS.SetBounds(0, N3PS.Size());
        N3PS[N3PS.Size()-1] = Agent1.p_sensor * Agent1.sensorweights[8];
        SS.SetBounds(0, SS.Size());
        SS[SS.Size()-1] = Agent1.state;
        N1SS.SetBounds(0, N1SS.Size());
        N1SS[N1SS.Size()-1] = Agent1.state * Agent1.sensorweights[3];
        N2SS.SetBounds(0, N2SS.Size());
        N2SS[N2SS.Size()-1] = Agent1.state * Agent1.sensorweights[6];
        N3SS.SetBounds(0, N3SS.Size());
        N3SS[N3SS.Size()-1] = Agent1.state * Agent1.sensorweights[9];

        Agent1.Step(BTStepSize, WorldFood);
        NO1.SetBounds(0, NO1.Size());
        NO1[NO1.Size()-1] = Agent1.NervousSystem.NeuronOutput(1);
        NO2.SetBounds(0, NO2.Size());
        NO2[NO2.Size()-1] = Agent1.NervousSystem.NeuronOutput(2);
        NO3.SetBounds(0, NO3.Size());
        NO3[NO3.Size()-1] = Agent1.NervousSystem.NeuronOutput(3);
        mov.SetBounds(0, mov.Size());
        mov[mov.Size()-1] = (Agent1.NervousSystem.NeuronOutput(2) - Agent1.NervousSystem.NeuronOutput(1));

        if (Agent1.birth == true){
            Agent1.state = Agent1.state - prey_repo;
            Agent1.birth = false;
        }
        if (Agent1.death == true){
            Agent1.Reset(rs.UniformRandomInteger(0,SpaceSize), 1.5);
            Agent1.death = false;
        }
        
        prey_pos.SetBounds(0, prey_pos.Size());
        prey_pos[prey_pos.Size()-1] = Agent1.pos;
        // Pred Sense & Step
        Agent2.Sense(prey_pos);
        Agent2.Step(BTStepSize, WorldFood, preylist);

        prey_pos.~TVector();
        pred_pos.~TVector();
        dead_food.~TVector();
    }
    // Save
    SSfile << FS << endl << N1FS << endl << N2FS << endl << N3FS << endl;
    SSfile << PS << endl << N1PS << endl << N2PS << endl << N3PS << endl;
    SSfile << SS << endl << N1SS << endl << N2SS << endl << N3SS << endl; 
    SSfile << NO1 << endl << NO2 << endl << NO3 << endl << mov << endl;

    SSfile.close();
}

// ---------------------------------------
// Eco Fear & Sensory Pollution Functions
// ---------------------------------------
void EcoFear(TVector<double> &genotype, RandomState &rs, double pred_condition)
{
    int maxpred = 30;
    int samplesize = 10;
    // Translate to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);
    // Create agents
    Prey prey(prey_netsize, prey_gain, prey_s_width, prey_frate, prey_feff, prey_metaloss, prey_movecost, prey_b_thresh);
    // Set nervous system
    prey.NervousSystem.SetCircuitSize(prey_netsize);
    int k = 1;
    // Prey Time-constants
    for (int i = 1; i <= prey_netsize; i++) {
        prey.NervousSystem.SetNeuronTimeConstant(i,phenotype(k));
        k++;
    }
    // Prey Biases
    for (int i = 1; i <= prey_netsize; i++) {
        prey.NervousSystem.SetNeuronBias(i,phenotype(k));
        k++;
    }
    // Prey Neural Weights
    for (int i = 1; i <= prey_netsize; i++) {
        for (int j = 1; j <= prey_netsize; j++) {
            prey.NervousSystem.SetConnectionWeight(i,j,phenotype(k));
            k++;
        }
    }
    // Prey Sensor Weights
    for (int i = 1; i <= prey_netsize*3; i++) {
        prey.sensorweights[i] = phenotype(k);
        k++;
    }
    ofstream EFfile("anallysis_results/EcoFear.dat");
    for (int j = 0; j<=maxpred; j++)
    {   
        printf("Collecting feeding rates with %d predators\n", j);
        TVector<double> lamb;
        for (int k = 0; k<=samplesize; k++){
            // Fill World to Carrying Capacity
            TVector<double> food_pos;
            TVector<double> WorldFood(1, SpaceSize);
            WorldFood.FillContents(0.0);
            for (int i = 0; i <= CC; i++){
                int f = rs.UniformRandomInteger(1,SpaceSize);
                WorldFood[f] = 1.0;
                food_pos.SetBounds(0, food_pos.Size());
                food_pos[food_pos.Size()-1] = f;
            }
            // Seed preylist with starting population
            TVector<Predator> predlist(0,-1);
            TVector<Prey> preylist(0,0);
            preylist[0] = prey;
            for (int i = 0; i < j; i++){
                Predator newpred(pred_gain, pred_s_width, pred_frate, pred_handling_time);
                newpred.Reset(rs.UniformRandomInteger(0,SpaceSize));
                newpred.frate = 0.0;
                newpred.condition = pred_condition;
                predlist.SetBounds(0, predlist.Size());
                predlist[predlist.Size()-1] = newpred;
                }

            double munch_count = 0;
            for (double time = 0; time < RateDuration; time += StepSize){
                // Remove chomped food from food list
                TVector<double> dead_food(0,-1);
                for (int i = 0; i < food_pos.Size(); i++){
                    if (WorldFood[food_pos[i]] <= 0){
                        dead_food.SetBounds(0, dead_food.Size());
                        dead_food[dead_food.Size()-1] = food_pos[i];
                    }
                }
                if (dead_food.Size() > 0){
                    for (int i = 0; i < dead_food.Size(); i++){
                        food_pos.RemoveFood(dead_food[i]);
                        food_pos.SetBounds(0, food_pos.Size()-2);
                    }
                }
                // Carrying capacity is 0 indexed, add 1 for true amount
                double food_count = food_pos.Size();
                double s_chance = G_Rate*food_count/(CC+1);
                double c = rs.UniformRandom(0,1);
                if (c > s_chance){
                    int f = rs.UniformRandomInteger(1,SpaceSize);
                    WorldFood[f] = 1.0;
                    food_pos.SetBounds(0, food_pos.Size());
                    food_pos[food_pos.Size()-1] = f;
                }
                // Prey Sense & Step
                TVector<double> pred_pos;
                for (int i = 0; i <= predlist.Size()-1; i++){
                    pred_pos.SetBounds(0, pred_pos.Size());
                    pred_pos[pred_pos.Size()-1] = predlist[i].pos;
                }
                prey.Sense(food_pos, pred_pos);
                prey.Step(StepSize, WorldFood);
                TVector<double> prey_pos(0,0);
                prey_pos[0] = prey.pos;
                // Check Births
                if (prey.birth == true){
                    prey.state = prey.state - prey_repo;
                    prey.birth = false;
                }
                // Check Deaths
                if (prey.death == true){
                    prey.Reset(rs.UniformRandomInteger(0,SpaceSize), 1.0);
                    prey.death = false;
                }
                // Check # of times food crossed
                munch_count += prey.snackflag;
                prey.snackflag = 0.0;  
                // Predator Sense & Step
                for(int i = 0; i < predlist.Size(); i++){
                    predlist[i].Sense(prey_pos);
                    predlist[i].Step(StepSize, WorldFood, preylist);
                }
                // Clear lists
                dead_food.~TVector();
                prey_pos.~TVector();
                pred_pos.~TVector();
            }

            double munchrate = munch_count/(RateDuration/StepSize);
            lamb.SetBounds(0, lamb.Size());
            lamb[lamb.Size()-1] = munchrate;
        }
        EFfile << lamb << endl;
        lamb.~TVector();
    }
    // Save
    EFfile.close();
}

void SPoll(TVector<double> &genotype, RandomState &rs, double pred_condition)
{
    // max * interval is greatest sensor width
    const int max_s = 30;
    const double s_interval = 10;
    const int samplesize = 10;
    // Translate to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);
    // Create agents
    Prey prey(prey_netsize, prey_gain, prey_s_width, prey_frate, prey_feff, prey_metaloss, prey_movecost, prey_b_thresh);
    // Set nervous system
    prey.NervousSystem.SetCircuitSize(prey_netsize);
    int k = 1;
    // Prey Time-constants
    for (int i = 1; i <= prey_netsize; i++) {
        prey.NervousSystem.SetNeuronTimeConstant(i,phenotype(k));
        k++;
    }
    // Prey Biases
    for (int i = 1; i <= prey_netsize; i++) {
        prey.NervousSystem.SetNeuronBias(i,phenotype(k));
        k++;
    }
    // Prey Neural Weights
    for (int i = 1; i <= prey_netsize; i++) {
        for (int j = 1; j <= prey_netsize; j++) {
            prey.NervousSystem.SetConnectionWeight(i,j,phenotype(k));
            k++;
        }
    }
    // Prey Sensor Weights
    for (int i = 1; i <= prey_netsize*3; i++) {
        prey.sensorweights[i] = phenotype(k);
        k++;
    }
    ofstream SPfile("analysis_results/SPoll.dat");
    for (int j = 0; j<=max_s; j++)
    {   
        prey.s_width = j*s_interval;
        printf("Collecting feeding rates with %f sensor width\n", prey.s_width);
        TVector<double> lamb;
        for (int k = 0; k<=samplesize; k++){
            // Fill World to Carrying Capacity
            TVector<double> food_pos;
            TVector<double> WorldFood(1, SpaceSize);
            WorldFood.FillContents(0.0);
            for (int i = 0; i <= CC; i++){
                int f = rs.UniformRandomInteger(1,SpaceSize);
                WorldFood[f] = 1.0;
                food_pos.SetBounds(0, food_pos.Size());
                food_pos[food_pos.Size()-1] = f;
            }
            // Seed preylist with starting population
            TVector<Predator> predlist(0,-1);
            TVector<Prey> preylist(0,0);
            preylist[0] = prey;
            for (int i = 0; i <= start_pred; i++){
                Predator newpred(pred_gain, pred_s_width, pred_frate, pred_handling_time);
                newpred.Reset(rs.UniformRandomInteger(0,SpaceSize));
                newpred.frate = 0.0;
                newpred.condition = pred_condition;
                predlist.SetBounds(0, predlist.Size());
                predlist[predlist.Size()-1] = newpred;
                }
            double munch_count = 0;
            for (double time = 0; time < RateDuration; time += StepSize){
                // Remove chomped food from food list
                TVector<double> dead_food(0,-1);
                for (int i = 0; i < food_pos.Size(); i++){
                    if (WorldFood[food_pos[i]] <= 0){
                        dead_food.SetBounds(0, dead_food.Size());
                        dead_food[dead_food.Size()-1] = food_pos[i];
                    }
                }
                if (dead_food.Size() > 0){
                    for (int i = 0; i < dead_food.Size(); i++){
                        food_pos.RemoveFood(dead_food[i]);
                        food_pos.SetBounds(0, food_pos.Size()-2);
                    }
                }
                // Carrying capacity is 0 indexed, add 1 for true amount
                double food_count = food_pos.Size();
                double s_chance = G_Rate*food_count/(CC+1);
                double c = rs.UniformRandom(0,1);
                if (c > s_chance){
                    int f = rs.UniformRandomInteger(1,SpaceSize);
                    WorldFood[f] = 1.0;
                    food_pos.SetBounds(0, food_pos.Size());
                    food_pos[food_pos.Size()-1] = f;
                }
                // Prey Sense & Step
                TVector<double> pred_pos;
                for (int i = 0; i <= predlist.Size()-1; i++){
                    pred_pos.SetBounds(0, pred_pos.Size());
                    pred_pos[pred_pos.Size()-1] = predlist[i].pos;
                }
                prey.Sense(food_pos, pred_pos);
                prey.Step(StepSize, WorldFood);
                TVector<double> prey_pos(0,0);
                prey_pos[0] = prey.pos;
                // Check Births
                if (prey.birth == true){
                    prey.state = prey.state - prey_repo;
                    prey.birth = false;
                }
                // Check Deaths
                if (prey.death == true){
                    prey.Reset(rs.UniformRandomInteger(0,SpaceSize), 1.0);
                    prey.death = false;
                }
                // Check # of times food crossed
                munch_count += prey.snackflag;
                prey.snackflag = 0.0;  
                // Predator Sense & Step
                for(int i = 0; i < predlist.Size(); i++){
                    predlist[i].Sense(prey_pos);
                    predlist[i].Step(StepSize, WorldFood, preylist);
                }
                // Clear lists
                dead_food.~TVector();
                prey_pos.~TVector();
                pred_pos.~TVector();
            }

            double munchrate = munch_count/(RateDuration/StepSize);
            lamb.SetBounds(0, lamb.Size());
            lamb[lamb.Size()-1] = munchrate;
        }
        SPfile << lamb << endl;
        lamb.~TVector();
    }
    // Save
    SPfile.close();
}

// ---------------------------------------
// Test Function for Code Development
// ---------------------------------------
void NewEco(TVector<double> &genotype, RandomState &rs)
{
    double test_CC = 1000;
    double test_frate = 0.5;
    double test_feff = 0.3;
    ofstream ppfile("analysis_results/TEST_prey_pop.dat");
    ofstream fpfile("analysis_results/TEST_food_pop.dat");
    // Translate to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);
    // Create agents
    // Playing with feff, frate, metaloss
    Prey Agent1(prey_netsize, prey_gain, prey_s_width, prey_frate, prey_feff, prey_metaloss, prey_movecost, prey_b_thresh);
    // Set nervous system
    Agent1.NervousSystem.SetCircuitSize(prey_netsize);
    int k = 1;
    // Prey Time-constants
    for (int i = 1; i <= prey_netsize; i++) {
        Agent1.NervousSystem.SetNeuronTimeConstant(i,phenotype(k));
        k++;
    }
    // Prey Biases
    for (int i = 1; i <= prey_netsize; i++) {
        Agent1.NervousSystem.SetNeuronBias(i,phenotype(k));
        k++;
    }
    // Prey Neural Weights
    for (int i = 1; i <= prey_netsize; i++) {
        for (int j = 1; j <= prey_netsize; j++) {
            Agent1.NervousSystem.SetConnectionWeight(i,j,phenotype(k));
            k++;
        }
    }
    // Prey Sensor Weights
    for (int i = 1; i <= prey_netsize*3; i++) {
        Agent1.sensorweights[i] = phenotype(k);
        k++;
    }
    // Fill World to Carrying Capacity
    TVector<double> food_pos;
    TVector<double> WorldFood(1, SpaceSize);
    WorldFood.FillContents(0.0);
    for (int i = 0; i <= CC; i++){
        int f = rs.UniformRandomInteger(1,SpaceSize);
        WorldFood[f] = 1.0;
        food_pos.SetBounds(0, food_pos.Size());
        food_pos[food_pos.Size()-1] = f;
    }
    // Make dummy predator list
    TVector<double> pred_pos(0,-1);
    TVector<Prey> preylist(0,0);
    preylist[0] = Agent1;
    // Carrying capacity is 0 indexed, add 1 for true amount
    for (int i = 0; i < 200; i++){
        double food_count = food_pos.Size();
        double s_chance = 1 - food_count/(test_CC+1);
        double c = rs.UniformRandom(0,1)*50;
        if (c < s_chance){
            int f = rs.UniformRandomInteger(1,SpaceSize);
            WorldFood[f] = 1.0;
            food_pos.SetBounds(0, food_pos.Size());
            food_pos[food_pos.Size()-1] = f;
        }
    }
    for (int i = 0; i < start_prey; i++){
        Prey newprey(prey_netsize, prey_gain, prey_s_width, test_frate, test_feff, prey_metaloss, prey_movecost, prey_b_thresh);
        newprey.Reset(rs.UniformRandomInteger(0,SpaceSize), 1.5);
        newprey.NervousSystem = Agent1.NervousSystem;
        newprey.sensorweights = Agent1.sensorweights;
        preylist.SetBounds(0, preylist.Size());
        preylist[preylist.Size()-1] = newprey;
        }
    for (double time = 0; time < PlotDuration; time += BTStepSize){
        // Remove chomped food from food list
        TVector<double> dead_food(0,-1);
        for (int i = 0; i < food_pos.Size(); i++){
            if (WorldFood[food_pos[i]] <= 0){
                dead_food.SetBounds(0, dead_food.Size());
                dead_food[dead_food.Size()-1] = food_pos[i];
            }
        }
        if (dead_food.Size() > 0){
            for (int i = 0; i < dead_food.Size(); i++){
                food_pos.RemoveFood(dead_food[i]);
                food_pos.SetBounds(0, food_pos.Size()-2);
            }
        }
        // Carrying capacity is 0 indexed, add 1 for true amount
        double food_count = food_pos.Size();
        double s_chance = 1 - food_count/(test_CC+1);
        double c = rs.UniformRandom(0,1)*50;
        if (c < s_chance){
            int f = rs.UniformRandomInteger(1,SpaceSize);
            WorldFood[f] = 1.0;
            food_pos.SetBounds(0, food_pos.Size());
            food_pos[food_pos.Size()-1] = f;
        }
        // Prey Sense & Step
        TVector<Prey> newpreylist;
        TVector<int> preydeaths;
        for (int i = 0; i < preylist.Size(); i++){
            preylist[i].Sense(food_pos, pred_pos);
            preylist[i].Step(BTStepSize, WorldFood);
            if (preylist[i].birth == true){
                preylist[i].state = preylist[i].state - prey_repo;
                preylist[i].birth = false;
                Prey newprey(prey_netsize, prey_gain, prey_s_width, test_frate, test_feff, prey_metaloss, prey_movecost, prey_b_thresh);
                newprey.NervousSystem = preylist[i].NervousSystem;
                newprey.sensorweights = preylist[i].sensorweights;
                newprey.Reset(preylist[i].pos+2, prey_repo);
                newpreylist.SetBounds(0, newpreylist.Size());
                newpreylist[newpreylist.Size()-1] = newprey;
            }
            if (preylist[i].death == true){
                preydeaths.SetBounds(0, preydeaths.Size());
                preydeaths[preydeaths.Size()-1] = i;
            }
        }
        // Update prey list with new prey list and deaths
        if (preydeaths.Size() > 0){
            for (int i = 0; i < preydeaths.Size(); i++){
                preylist.RemoveItem(preydeaths[i]);
                preylist.SetBounds(0, preylist.Size()-2);
            }
        }
        if (newpreylist.Size() > 0){
            for (int i = 0; i < newpreylist.Size(); i++){
                preylist.SetBounds(0, preylist.Size());
                preylist[preylist.Size()-1] = newpreylist[i];
            }
        }
        ppfile << preylist.Size() << endl;
        fpfile << food_pos.Size() << endl;
        // Check Population Collapse
        if (preylist.Size() <= 0){
            break;
        }
        else{
            newpreylist.~TVector();
            preydeaths.~TVector();
            dead_food.~TVector();
        }
    }
    // Save
    ppfile.close();
    fpfile.close();
}

// ================================================
// E. MAIN FUNCTION
// ================================================
int main (int argc, const char* argv[]) 
{
// ================================================
// EVOLUTION
// ================================================
	// TSearch s(VectSize);
    // long seed = static_cast<long>(time(NULL));
	// // save the seed to a file
	// ofstream seedfile;
	// seedfile.open ("seed.dat");
	// seedfile << seed << endl;
	// seedfile.close();
	// // Configure the search
	// s.SetRandomSeed(seed);
	// s.SetSearchResultsDisplayFunction(ResultsDisplay);
	// s.SetPopulationStatisticsDisplayFunction(EvolutionaryRunDisplay);
	// s.SetSelectionMode(RANK_BASED);
	// s.SetReproductionMode(GENETIC_ALGORITHM);
	// s.SetPopulationSize(POPSIZE);
	// s.SetMaxGenerations(GENS);
	// s.SetCrossoverProbability(CROSSPROB);
	// s.SetCrossoverMode(UNIFORM);
	// s.SetMutationVariance(MUTVAR);
	// s.SetMaxExpectedOffspring(EXPECTED);
	// s.SetElitistFraction(ELITISM);
	// // s.SetSearchConstraint(1);
    // ofstream evolfile;
	// evolfile.open("fitness.dat");
	// cout.rdbuf(evolfile.rdbuf());

    // // Staged Evolution: Progressively Reduce Carrying Capacity
    // for (int i = maxCC; i >= minCC; i--){
    //     int CC = i;
    //     printf("Carrying Capacity: %d\n", CC);
    //     s.SetSearchTerminationFunction(CCTerminationFunction);
    //     s.SetEvaluationFunction(Coexist);
    //     s.ExecuteSearch();
    //     ofstream S1B;
    //     TVector<double> bestVector1 = s.BestIndividual();
    //     S1B.open("best.gen.dat");
    //     cout.rdbuf(S1B.rdbuf());
    //     S1B << setprecision(32);
    //     S1B << bestVector1 << endl;
    //     S1B.close();
    //     cout.rdbuf(evolfile.rdbuf());
    // }

    // printf("Final Tuning\n");
    // int CC = finalCC;
    // s.SetSearchTerminationFunction(EndTerminationFunction);
    // s.SetEvaluationFunction(Coexist);
    // s.ExecuteSearch();
    // ofstream S1B;
    // TVector<double> bestVector1 = s.BestIndividual();
    // S1B.open("best.gen.dat");
    // cout.rdbuf(S1B.rdbuf());
    // S1B << setprecision(32);
    // S1B << bestVector1 << endl;
    // S1B.close();
    // cout.rdbuf(evolfile.rdbuf());
    
    // return 0;

// ================================================
// RUN ANALYSES
// ================================================
    // load the seed
    ifstream seedfile;
    double seed;
    seedfile.open("seed.dat");
    seedfile >> seed;
    seedfile.close();
    // load best individual
    ifstream genefile;
    genefile.open("menagerie/bigredgen2.dat");
    TVector<double> genotype(1, VectSize);
	genefile >> genotype;
    genefile.close();
    // set the seed
    RandomState rs(seed);

    // // Behavioral Traces // // 
    // BehavioralTracesCoexist(genotype, rs, pred_condition);

    // // Interaction Rate Collection // //
    // CollectEcoRates(genotype, rs);

    // // Sensory Sample Collection // //
    // SensorySample(genotype, rs);

    // // EcoFear Analysis // // 
    // EcoFear(genotype, rs, pred_condition);

    // // Sensory Pollution Analysis // // 
    // SPoll(genotype, rs, pred_condition);

    // // Code Testbed // // 
    NewEco(genotype, rs);

    return 0;

}