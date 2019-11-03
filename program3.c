// Daniel Patrick
// CS 620 Programming Assignment 3 - Spring 2018

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Track how often the compute_fitness function is called
int fitness_count=0;


/* **** STRUCT **** */
// data_item_struct provide weight and profit for an item of knapsack problem
struct data_item_struct {
  int weight;
  int profit;
};
typedef struct data_item_struct data_item;

// data_knapsack_struct provide capacity, count of item, and set of items of knapsack problem
struct data_knapsack_struct {
  int capacity;
  int size;
  data_item* item;
};
typedef struct data_knapsack_struct data_knapsack;

// chromosome_struct provide set of genes of a chromosome and fitness value of the chromosome
struct chromosome_struct{
  int* gene;
  int fitness;
};
typedef struct chromosome_struct chromosome;


/* **** DATA FUNCTION **** */
// Random real number between min and max values
double random_range(double min, double max);

// Read knapsack data problem from <filename>, put them into <data>, which must already be allocated
void read_data(char* filename, data_knapsack* data);

// Find the best fitness value in generation <gen>
int find_best_fitness(chromosome* pop, int pop_size, int gen);

// Print weight and profit of the <item>
void print_item(data_item * item);

// Print weight and profit of all items
void print_items(data_knapsack data);

// Print capacity, count of item, and all items from knapsack problem
void print_data(data_knapsack data);

// Print gene and fitness value of the chromosome<pop>
void print_chromosome(chromosome pop, int gene_size);

// Print gene and fitness value of all chromosomes in the generation<gen>
void print_chromosomes(chromosome* pop, int gene_size,int pop_size, int gen);

// Print gene and fitness value of all chromosomes in every generation
void print_all_chromosomes(chromosome* pop, int gene_size,int pop_size, int max_gen);

// Allocate gene that has length <gene_size> for the chromosome <c>
void allocate_gene(int gene_size, chromosome* c);

// Allocate gene that has length <gene_size> for all chromosomes in every generation
void allocate_all_genes(chromosome* pop, int gene_size, int pop_size, int gen_size);

// Allocate all chromosomes
// Return pointer to <population>
chromosome* allocate_all_chromosomes(int gene_size, int pop_size, int gen_size);

// Random turn off one gene from the chromosome that is overweight
void turn_off_one_bit(int* solution, int count);

// Compute total fitness for all chromosomes in the generation<gen>
int compute_total_fitness(chromosome* pop, int pop_size, int gen);

// Compute and return fitness value (total profit of <gene>)
// Fitness value is the total value of profit of all gene that has value 1
int compute_fitness(data_knapsack* data, int* gene);

// Initialize the chromosome <c>
// Random 0 or 1 for each gene in the chromosome <c>
// Compute fitness of the chromosome <c>
void initialize_random_chromosome(data_knapsack* data, chromosome* c);

// Initialize all the chromosomes for the first generation
void initialize_population(data_knapsack* data, int pop_size, chromosome* pop);

// ROULETE WHEEL SELECTION
// Random select the chromosome by using roulette wheel
// Roulette wheel is random choose the chromosome by the percentage of chromosome fitness
// Return the position that is selected
int roulette_wheel_selection(chromosome* pop, int pop_size, int total_fitness, int gen);

// CROSSOVER
// Crossover chromosomes from positions <from_a> and <from_b>, 
// and put the new chromosomes to positions <to_a> and <to_b>
// Compute the fitness value of the new chromosomes
void crossover_chromosomes(chromosome *pop, int from_a, int from_b, int to_a, int to_b, int gene_size, data_knapsack* data);

// MUTATION
// Mutate chromosome in the position <which>
// Random switch genes by the percentage <mutation_pct>
// Compute the fitness value if it has mutated
void mutate_chromosome(chromosome *pop, int which, int gene_size, double mutation_pct, data_knapsack* data);

// COPY
// Copy chromosome <from> to chromosome <to>
void copy_chromosomes(chromosome *pop, int from, int to, int gene_size);

// Set variables that refer to parent position and child position
// Use roulette wheel selection for finding the parent position
void set_variables(chromosome* pop, int pop_size, int total_fitness, int gen, int* from_a, int* from_b, int* to_a, int* to_b, int i);

/* Include the implementation of the functions that are defined above */
#include "program3-lib.c"


/* **** MAIN FUNCTION **** */
int main(int argc, char ** argv) {
  // Set up the seed random for random function
  srandom(time(NULL));

  // Verify that the correct number of inputs were provided
  if ( argc < 6 ) {
    printf("Usage: %s FILENAME POPULATION_SIZE CROSSOVER_PCT MUTATION_PCT MAX_GENERATIONS\n", argv[0]);
    exit(1);
  }

  // Set the percentage of crossover and mutation
  double crossover_pct = atof(argv[3]);
  double mutation_pct = atof(argv[4]);
  // Set the population size and the number of generation
  int population_size = atoi(argv[2]);
  int max_generation = atoi(argv[5]);
  
  // Print out the command line for verification
  printf("Filename: %s\tPopulation Size: %d\tCrossover PCT: %f\tMutation PCT: %f\tMax Genaration: %d\n", argv[1], population_size, crossover_pct, mutation_pct, max_generation);

  // Allocate the initial data for knapsack problem
  data_knapsack data;

  // Read knapsack problem data from file <Filename>
  read_data(argv[1], &data);
  // Print the knapsack problem data
  print_data(data);

  // Allocate the population for all chromosomes in every generation
  chromosome* population = allocate_all_chromosomes(data.size, population_size, max_generation);

  // Initialize chromosomes in the first generation
  initialize_population(&data, population_size, population);

  int gen;// Current generation
  // Loop each generation untill the total number of generation - 1
  // because the parent of the last generation compute the last generation
  for(gen=0 ; gen < max_generation-1 ; gen++){
    // Set total fitness that compute from all chromosomes in the current generation
    // Total fitness is used for roulette wheel selection
    int total_fitness = compute_total_fitness(population, population_size, gen);

    int i;
    for(i=0; i < population_size; i+=2){
      int from_a, from_b, to_a, to_b;
      // Set variables that refer to parent position <from_a, from_b> and child position <to_a, to_b>
      // Use roulette wheel selection for finding the parent position
      set_variables(population, population_size, total_fitness, gen, &from_a, &from_b, &to_a, &to_b, i);
      // Do crossover by the percentage <crossover_pct>
      if(random_range(0.0, 1.0) < crossover_pct){
        crossover_chromosomes(population, from_a, from_b, to_a, to_b, data.size, &data);
      }
      else{
        copy_chromosomes(population, from_a, to_a, data.size);
        copy_chromosomes(population, from_b, to_b, data.size);
      }
      mutate_chromosome(population, to_a, data.size, mutation_pct, &data);
      mutate_chromosome(population, to_b, data.size, mutation_pct, &data);
    }
  }

  print_all_chromosomes(population, data.size, population_size, max_generation);
  printf("\nBest fitness: %d\n", population[find_best_fitness(population, population_size, max_generation-1)].fitness);
  printf("Fitness count: %d\n", fitness_count);
  exit(0);
}
