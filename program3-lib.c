// Daniel Patrick
// CS 620 Programming Assignment 3 - Spring 2018

// Random real number between min and max values
double random_range(double min, double max){
  return min + (random() / ( RAND_MAX / (max-min) ) ) ;
}

// Read knapsack data problem from <filename>, put them into <data>, which must already be allocated
void read_data(char * filename, data_knapsack* data) {

  FILE * input_file;

  input_file = fopen(filename, "r");

  if ( input_file == NULL ) {
    printf("Cannot open input file: %s\n", filename);
    exit(1);
  }

  fscanf(input_file, "%d", &(data->capacity));
  fscanf(input_file, "%d", &(data->size));

  data->item = malloc(data->size * sizeof(data_item));

  int n = 0;
  int retval;

  while ( n < data->size ) {
    if ( feof(input_file) ) {
      printf("Not enough lines in data file\n");
      exit(1);
    }

    retval = fscanf(input_file, "%d %d",
        &(data->item[n].weight),
        &(data->item[n].profit));

    if ( retval < 2 ) {
      printf("Not enough lines in data file\n");
      exit(1);
    }

    n++;
  }
}

// Find the best fitness value in generation <gen>
int find_best_fitness(chromosome* pop, int pop_size, int gen){
  int max = gen * pop_size;
  int i;
  for(i=0 ; i < pop_size ; i++){
    if(pop[max].fitness < pop[gen * pop_size + i].fitness){
      max = gen * pop_size + i;
    }
  }
  return max;
}


/* **** PRINT FUNCTION **** */
// Print weight and profit of the <item>
void print_item(data_item * item) {
  printf("Weight: %d\t Profit: %d\n", item->weight, item->profit);
}
// Print weight and profit of all items
void print_items(data_knapsack data) {
  int i;
  for ( i = 0; i < data.size; i++) {
    print_item(&data.item[i]);
  }
}
// Print capacity, count of item, and all items from knapsack problem
void print_data(data_knapsack data) {
  printf("Capacity: %d\n", data.capacity);
  printf("Size: %d\n", data.size);

  printf("--- Data:\n");
  print_items(data);
}

// Print gene and fitness value of the chromosome<pop>
void print_chromosome(chromosome pop, int gene_size){
  int i;
  for(i=0 ; i < gene_size ; i++){
    printf("%d ", pop.gene[i]);
  }
  printf("\tFitness: %d\n", pop.fitness);
}
// // Print gene and fitness value of all chromosomes in the generation<gen>
void print_chromosomes(chromosome* pop, int gene_size,int pop_size, int gen){
  int i;
  for(i=0 ; i < pop_size ; i++){
    printf("Chromosome %d: ", i);
    print_chromosome(pop[gen * pop_size + i], gene_size);
  }
}
// Print gene and fitness value of all chromosomes in every generation
void print_all_chromosomes(chromosome* pop, int gene_size,int pop_size, int max_gen){
  int i;
  for(i=0 ; i < max_gen ; i++){
    printf("\nGenaration %d: \n", i);
    print_chromosomes(pop, gene_size, pop_size, i);
  }
}


/* **** ALLOCATE FUNCTION **** */
// Allocate gene that has length <gene_size> in the chromosome <c>
void allocate_gene(int gene_size, chromosome* c){
  c->gene = malloc(gene_size * sizeof(int));
}
// Allocate gene that has length <gene_size> for all chromosomes in every generation
void allocate_all_genes(chromosome* pop, int gene_size, int pop_size, int gen_size){
  int i;
  for(i=0 ; i < pop_size * gen_size; i++){
    allocate_gene(gene_size, &pop[i]);
  }
}
// Allocate all chromosomes
// Return pointer to <population>
chromosome* allocate_all_chromosomes(int gene_size, int pop_size, int gen_size){
  chromosome* population = malloc(pop_size * gen_size * sizeof(chromosome));
  allocate_all_genes(population, gene_size, pop_size, gen_size);
  return population;
}


/* **** GENETIC FUNCTION **** */
// Random turn off one gene from the chromosome that is overweight
void turn_off_one_bit(int* solution, int count){
  int c = 0;
  int where[count];
  int i;

  // Find the positions that have 1 not 0
  for(i=0 ; i < count ; i++){
    if(solution[i] == 1){
      where[c] = i;
      c++;
    }
  }

  int j = random() % c; // Random only the positions that have 1
  solution[where[j]] = 0; // Set the random position to 0
}
// Compute total fitness from all chromosomes in the generation(gen)
int compute_total_fitness(chromosome* pop, int pop_size, int gen){
  int total_fitness = 0;
  int i;
  for(i=0 ; i < pop_size ; i++){
    total_fitness += pop[gen * pop_size + i].fitness;
  }
  return total_fitness;
}
// Compute and return fitness value (total profit of <gene>)
// Fitness value is the total value of profit of all gene that has value 1
int compute_fitness(data_knapsack* data, int* gene){
  int weight=0;
  int profit=0;
  fitness_count++; // Track how often the compute_fitness function is called
  int i;
  for(i=0 ; i < data->size ; i++){
    if(gene[i] == 1){
      weight += data->item[i].weight;
      profit += data->item[i].profit;
    }
  }
  // If total weight of the chromosome more than capacity, 
  // random turn of one bit of gene and compute fitness again
  if(weight > data->capacity){
    turn_off_one_bit(gene, data->size);
    return compute_fitness(data, gene); // Call compute_fitness function recursively
  }

  // Return the total profit (fitness value) of all gene <gene>
  return profit;
}
// Initialize the chromosome <c>
// Random 0 or 1 for each gene in the chromosome <c>
// Compute fitness of the chromosome <c>
void initialize_random_chromosome(data_knapsack* data, chromosome* c){
  int i;
  // Random 0 or 1 for each gene in the chromosome
  for(i=0 ; i < data->size ; i++){
    c->gene[i] = random() % 2;
  }
  // Compute fitness of the chromosome <c>
  c->fitness = compute_fitness(data, c->gene);
}
// Initialize all the chromosomes for the first generation
void initialize_population(data_knapsack* data, int pop_size, chromosome* pop){
  int i;
  for(i=0 ; i < pop_size; i++){
    initialize_random_chromosome(data, &pop[i]);
  }
}
// ROULETE WHEEL SELECTION
// Random select the chromosome by using roulette wheel
// Roulette wheel is random choose the chromosome by the percentage of chromosome fitness
// Return the position that is selected
int roulette_wheel_selection(chromosome* pop, int pop_size, int total_fitness, int gen){
  int chunk = random() % total_fitness; // Random number less than <total_fitness>
  int i;
  for(i=0 ; i < pop_size; i++){
    if(chunk <= pop[gen*pop_size + i].fitness){
      return i; // Return the position that is selected
    }
    chunk = chunk - pop[gen*pop_size + i].fitness;
  }

  exit(1); // should never reach this point
  return 0;
}
// CROSSOVER
// Crossover chromosomes from positions <from_a> and <from_b>, 
// and put the new chromosomes to positions <to_a> and <to_b>
// Compute the fitness value of the new chromosomes
void crossover_chromosomes(chromosome *pop, int from_a, int from_b, int to_a, int to_b, int gene_size, data_knapsack* data){
  int spot = random() % (gene_size - 1); // Random number for selecting the position to crossover
  int i;
  // Copy genes from <from_a> to <to_a> and <from_b> to <to_b> until the position "spot"
  for (i=0 ; i <= spot ; i++){
    pop[to_a].gene[i] = pop[from_a].gene[i];
    pop[to_b].gene[i] = pop[from_b].gene[i];
  }
  // Copy genes from position "spot" to the end from <from_b> to <to_a> and <from_a> to <to_b>
  for (i=spot+1 ; i < gene_size; i++){
    pop[to_a].gene[i] = pop[from_b].gene[i];
    pop[to_b].gene[i] = pop[from_a].gene[i];
  }
  // Compute the fitness value of the new chromosomes
  pop[to_a].fitness = compute_fitness(data, pop[to_a].gene);
  pop[to_b].fitness = compute_fitness(data, pop[to_b].gene);
}
// MUTATION
// Mutate chromosome in the position <which>
// Random switch genes by the percentage <mutation_pct>
// Compute the fitness value if it has mutated
void mutate_chromosome(chromosome *pop, int which, int gene_size, double mutation_pct, data_knapsack* data){
  int i;
  int c = 0;
  for (i=0 ; i < gene_size ; i++){
    // Random switch genes by the percentage <mutation_pct>
    if ( random_range(0.0, 1.0) < mutation_pct){
      pop[which].gene[i] = 1 - pop[which].gene[i];
      c++;
    }
  }
  // Compute the fitness value if it has mutated
  if (c > 0){
    pop[which].fitness = compute_fitness(data, pop[which].gene);
  }

}
// COPY
// Copy chromosome <from> to chromosome <to>
void copy_chromosomes(chromosome *pop, int from, int to, int gene_size){
  int i;
  // Copy genes
  for (i=0 ; i < gene_size ; i++){
    pop[to].gene[i] = pop[from].gene[i];
  }
  // Copy fitness
  pop[to].fitness = pop[from].fitness;
}
// Set variables that refer to parent position and child position
// Use roulette wheel selection for finding the parent position
void set_variables(chromosome* pop, int pop_size, int total_fitness, int gen, int* from_a, int* from_b, int* to_a, int* to_b, int i){
  *from_a = gen*pop_size + roulette_wheel_selection(pop, pop_size, total_fitness, gen);
  *from_b = gen*pop_size + roulette_wheel_selection(pop, pop_size, total_fitness, gen);
  *to_a = ((gen+1) * pop_size) + i;
  *to_b = ((gen+1) * pop_size) + i+1;

  // Check <to_be> is refer to position that is not in the same generation
  // It will be happen when the population size is odd number
  if(*to_b >= (gen+2) * pop_size){
    *to_b = *to_a;
  }
}