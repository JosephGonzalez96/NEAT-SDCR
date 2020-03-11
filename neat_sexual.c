#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <stdbool.h>
#include <string.h>
// #include <omp.h>
// #include <mpi.h>
#include <assert.h>


// int population_size = 100;
// float population_mortality = 0.125;
float population_mortality = 0.75;
float distance_compatibility_in_specie = 30;
// float compatible_fitness = 3;
//----------------------------------------------------------
float chance_mutation_add_node = 0.2;
float chance_mutation_add_connection = 0.3;
float chance_mutation_change_weights = 0.9;
//----------------------------------------------------------
float chance_mutation_weight_per_genome = 0.7;
float variation_weights = 0.1;
float variation_bias = 0.1;
float variation_activation = 0.1;
//------------------------------------------------------------
int age_importance = 2;
int advance_in_fitness = 0;
int generation_stop = 1200;
float c1 = 0.1, c2 = 4.0, N=1;
int tests_numbers = 100;

//-------------------------------------------------------------------------
float structural_connections[2] = {0.6,0.90}; //aqui va a acumulado //hid-hid, ext-hid, ext-ext
// cambios en la funcion: calculate offspring, select partner, sexual reproduction
// no tocar estos parametros
int damned;
int innovation = 1;
int kills;
double total_fitness;
int number_species = 1;
int living_species = 0;
int epsilon = 0;
int new_childs;
int reproductive_species;
// variables del auto 
int max_velocity = 10;              // en cm/s
float reaction_time = 1;            // en segundos
int car_length = 10;                // en centimetros
int max_angle = 45;                 // en grados 
// variables de la imagen
//225 mal llega solo hasta 224 con h=96, l=128
int number_input_images = 250;
// int height = 480;                //cambiar abajo tambien el tamaño
// int length = 640;
// int height = 96;                //cambiar abajo tambien el tamaño
// int length = 128;
int height = 2;             
int length = 2;
int number_inputs = 2*2*4;
float color_propability = 0.8;
int number_outputs = 2;
// population 
int max_generations = 300;
int population = 500;
int number_offspring = 100;
int zeros;

int maxim_name = 1;

clock_t start, end;
double cpu_time_used;
int average_size;

// --------------------------------Estructuras----------------------------------
typedef struct Gene{
    int into;
    int out;
    int origin;
    float weight; 
    struct Gene *next;
    bool is_bias;
    float bias;
    bool is_activation;
    float activation;
} Gene;
typedef struct Genome{
    struct Gene *head;
    struct Neuron *neurons;
    struct Neuron *last_neuron;
    struct Genome *next;
    struct Genome *previous;
    int hiden_neurons;
    int name;
    int size;
    int sons;
    int age;
    int rank;
    double fitness;
    bool death_sentence;
    bool occupied_outputs[2];
} Genome;
typedef struct Eval_tree{
    int neuron;
    double value;
    float activation;
    struct Eval_tree *left_tree;
    struct Eval_tree *right_tree;
    struct Eval_tree *parent_tree;
} Eval_tree; 
typedef struct Specie{
    int identifier;
    struct Genome *first;
    struct Genome *last;
    struct Specie *next;
    struct Distances *distances;
    double fitness;
    double max_fitness;
    double previous_max_fitness;
    int it_advance;
    int population;
    int reproductions;
    int new_childs;
    int max_name;
    int reproductivity;
} Specie;
typedef struct Distances{
    int identifier;
    float distance;
    struct Distances *next;
} Distances;
typedef struct Neuron{
    int identifier;
    struct Neuron *next;
} Neuron;
//----------------------------------Funciones para cargar la imagen-------------------------------
void charge_inputs_and_outputs(unsigned char (*inputs)[height][length][3],int (*outputs)[length]){
    char name_file_image[25];
    char name_file_points[25];
    char number[5];
    char line[15*length];
    char *c_line;
    char *phrase;
    char *data;
    int max_number;
    FILE *file_image, *file_points;
    for(int k=0; k<number_input_images; k++ ){
        strcpy(name_file_image,"Pixels/pixels_");
        sprintf(number,"%d",k+1); 
        strcat(name_file_image,number);
        strcat(name_file_image,".txt");
        file_image = fopen(name_file_image,"r");
        for(int i=0; i<height; i++){
            fgets (line ,15*length,file_image);                          //para obtener cada linea
            c_line = line;
            for(int j=0; j<length; j++){
                phrase = strtok_r(c_line,";",&c_line);
                inputs[k][i][j][0] = (unsigned char) atoi(strtok_r(phrase,",",&phrase));
                inputs[k][i][j][1] = (unsigned char) atoi(strtok_r(phrase,",",&phrase));
                inputs[k][i][j][2] = (unsigned char) atoi(strtok_r(phrase,",",&phrase));
            }
        }
        fclose(file_image);
        strcpy(name_file_points,"Points/points_");
        strcat(name_file_points,number);
        strcat(name_file_points,".txt");
        file_points = fopen(name_file_points,"r");
        fgets (line ,15*length,file_points);
        c_line = line;
        max_number = atoi(strtok_r(c_line,";",&c_line));
        if(max_number == 0){
            max_number = 1;
        }
        for(int j=0; j<length; j++){
            outputs[k][j] = ((double) atoi(strtok_r(c_line,",",&c_line))/(double) max_number)*255;
        }  
        fclose(file_points);
        printf("cargado: %d\n",k+1);
    }
}
// ---------------------------------Funciones--------------------------------------
void find_coordenates(int number, int (*coordinates)){
    int copy_number = number;
    int i = 0, j = 0;
    int coor[3];
    while(copy_number > length*4){
        copy_number += -length*4;
        i++;
    }
    while(copy_number > 4){
        copy_number += -4;
        j++;
    }  
    coordinates[0]=i;
    coordinates[1]=j;
    coordinates[2]=copy_number;
}
void print_all_neurons(Genome *genome){
    struct Neuron *current_neuron = genome->neurons;
    if(genome->neurons!=NULL){
        printf("hidden neurons:%d, first:%d, last:%d\n",genome->hiden_neurons,genome->neurons->identifier,genome->last_neuron->identifier);
    }else{
        printf("hidden neurons:%d\n",genome->hiden_neurons);
    }
    for(int i=0; i<genome->hiden_neurons; i++){
        printf("id_neuron:%d \n",current_neuron->identifier);
        current_neuron = current_neuron->next;
    }
}
void print_genome(Genome* genome){
    struct Gene* gene = genome->head;
    printf("The name is %d and the fitness is: %f and age: %d and reproductivity: %d \n",genome->name,genome->fitness,genome->age,genome->rank + genome->age);
}
void print_genome_genes(Genome* genome, int size){
    struct Gene* gene = genome->head;
    int coordinates[3];
    printf("\nThe name is: %d, hidden neurons:%d, size:%d\n",genome->name,genome->hiden_neurons,genome->size);
    for(int i=1; i<=genome->size; i++){
        // printf("Gene %d -> into:%d, outo:%d, weight:%f, origin:%d, bias:%f, is_bias:%d\n",i,gene->into,gene->out,gene->weight,gene->origin,gene->bias,gene->is_bias)
        printf("Gene %d -> ",i);
        if(gene->into <= number_inputs){
            find_coordenates(gene->into,coordinates);
            printf("into:(%d,%d,%d),",coordinates[0]+1,coordinates[1]+1,coordinates[2]);
        }else{
            printf("into:%d,",gene->into);
        }
        if(gene->out>number_inputs && gene->out<=number_inputs+number_outputs){
            printf(" outo:%d, ",gene->out-number_inputs);
        }else{
            printf(" outo:%d, ",gene->out);
        }
        printf("weight:%f, origin:%d, bias:%f, is_bias:%d\n",gene->weight,gene->origin,gene->bias,gene->is_bias);
        gene = gene->next;
    }
}
void print_genome_in_files(Genome* genome, FILE *file_genomes){
    struct Gene* gene = genome->head;
    int coordinates[3];
    fprintf(file_genomes,"\nThe name is: %d, hidden neurons:%d, size:%d\n",genome->name,genome->hiden_neurons,genome->size);
    for(int i=1; i<=genome->size; i++){
        fprintf(file_genomes,"Gene %d -> into:%d, outo:%d, weight:%.16f, origin:%d, bias:%.16f, is_bias:%d, acti:%.16f, is_acti:%d \n",i,gene->into,gene->out,gene->weight,gene->origin,gene->bias,gene->is_bias,gene->activation,gene->is_activation);
        
        // fprintf(file_genomes,"Gene %d -> ",i);
        // if(gene->into <= number_inputs){
        //     find_coordenates(gene->into,coordinates);
        //     fprintf(file_genomes,"into:(%d,%d,%d),",coordinates[0]+1,coordinates[1]+1,coordinates[2]);
        // }else{
        //     fprintf(file_genomes,"into:%d,",gene->into);
        // }
        // if(gene->out>number_inputs && gene->out<=number_inputs+number_outputs){
        //     fprintf(file_genomes," outo:%d, ",gene->out-number_inputs);
        // }else{
        //     fprintf(file_genomes," outo:%d, ",gene->out);
        // }
        // fprintf(file_genomes,"weight:%.15f, origin:%d, bias:%.15f, is_bias:%d\n",gene->weight,gene->origin,gene->bias,gene->is_bias);
        
        gene = gene->next;
    }
    
    // if(gene->into <= number_inputs){
    //     find_coordenates(gene->into,coordinates);
    //     fprintf(file_evaluation,"Gene:%d -> in:(%d,%d,%d)%d  out:%d\n",i,coordinates[0]+1,coordinates[1]+1,coordinates[2],gene->into,gene->out);
    // }else{
    //     fprintf(file_evaluation,"Gene:%d -> in:%d  out:%d\n",i,gene->into,gene->out);
    // }
}
void print_species(Specie *specie){
    struct Genome *first_genome;
    struct Genome *current_genome;
    while(specie!=NULL){
        first_genome = specie->first;
        current_genome = first_genome;
        for(int i=1; i<=specie->population; i++){
            printf("The name is %d and the fitness is: %d and age: %d and reproductivity: %d \n",current_genome->name,(int)current_genome->fitness,current_genome->age,current_genome->rank + current_genome->age);
            // print_genome(current_genome);
            print_genome_genes(current_genome,current_genome->size);
            current_genome = current_genome->next;
        }
        // printf("----------------------------------------------------------------------------------\n");
        specie = specie->next;
    }
}
//-----------------------------------Funciones para las neuronas------------------------------------------------
Neuron* create_new_neuron(int identifier){
    struct Neuron *neuron = malloc(sizeof(Neuron));
    neuron->identifier = identifier;
    neuron->next = NULL;
    return neuron;
}
void hitch_neuron_to_genome(int neuron_identifier, Genome *genome){
    struct Neuron *neuron = create_new_neuron(neuron_identifier);
    if(genome->last_neuron != NULL){
        genome->last_neuron->next = neuron;
    }else{
        genome->neurons = neuron;
    }
    genome->last_neuron = neuron;
}
bool is_neuron_present_in_genome(int identifier_neuron, Genome *genome){
    struct Neuron *current_neuron = genome->neurons;
    bool finded = false;
    while(!finded && current_neuron!=NULL){
        if(current_neuron->identifier == identifier_neuron){
            finded = true;
        }
        current_neuron = current_neuron->next;
    }
    return finded;
}
void copy_all_neurons(Genome *parent, Genome *child){
    struct Neuron *neuron = parent->neurons;
    for(int i=0; i<parent->hiden_neurons; i++){
        hitch_neuron_to_genome(neuron->identifier,child);
        neuron = neuron->next;
    }
}
void copy_all_neurons_for_mix_genomes(Genome *parent1,Genome *parent2, Genome *child){
    struct Neuron *neuron = parent1->neurons;
    int hidden = 0;
    for(int i=0; i<parent1->hiden_neurons; i++){
        hitch_neuron_to_genome(neuron->identifier,child);
        hidden++;
        neuron = neuron->next;
    }
    neuron = parent2->neurons;
    for(int i=0; i<parent2->hiden_neurons; i++){
        if(!is_neuron_present_in_genome(neuron->identifier,child)){
            hitch_neuron_to_genome(neuron->identifier,child);
            hidden++;
        }
        neuron = neuron->next;
    }
    child->hiden_neurons = hidden;
}
void delete_all_neurons(Genome *genome){
    struct Neuron *neuron = genome->neurons;
    struct Neuron *next;
    for(int i=0; i<genome->hiden_neurons; i++){
        next = neuron->next;
        free(neuron);
        neuron = next;
    }
}
//----------------------------------Funciones para la mutacion y manipulacion de genes-----------------------
Genome* create_new_genome(){
    struct Genome *genome =  malloc(sizeof (Genome));
    genome->head = NULL;
    genome->next = NULL;
    genome->previous = NULL;
    genome->neurons = NULL;
    genome->last_neuron = NULL;
    genome->hiden_neurons = 0;
    genome->name = maxim_name;
    genome->size = 0;
    genome->sons = 0;
    genome->age = 0;
    genome->rank = 0;
    genome->fitness = 0;
    genome->death_sentence = false;
    genome->occupied_outputs[0]=false;
    genome->occupied_outputs[1]=false;
    return genome; 
}
Gene* create_new_gene(int into, int out, int weight_code, float weight, int bias_code,float bias,bool is_bias,int activation_code,float activation, bool is_activation, int origin, struct Gene* next){
    struct Gene *new_gene = malloc(sizeof(Gene));
    new_gene->into = into; 
    new_gene->out = out;
    if(weight_code == 1){
        new_gene->weight = weight;
    }else if(weight_code == 2){
        new_gene->weight = 1;
    }else if(weight_code == 3){
        new_gene->weight = (((float)rand()/(float)RAND_MAX)*2)-1;
    }
    if(bias_code == 1){
        new_gene->bias = bias;
        new_gene->is_bias = is_bias;
    }else if(bias_code == 2){
        new_gene->bias = 0;
        new_gene->is_bias = false;
    }else if(bias_code == 3){
        new_gene->bias = (((float)rand()/(float)RAND_MAX)*2)-1;
        new_gene->is_bias = true;
    }
    if(activation_code == 1){
        new_gene->activation = activation;
        new_gene->is_activation = is_activation;
    }else if(activation_code == 2){
        new_gene->activation = 0;
        new_gene->is_activation = false;
    }else if(activation_code == 3){
        new_gene->activation = 0.5;
        new_gene->is_activation = true;
    }
    new_gene->origin = origin;
    new_gene->next = next;
    return new_gene;
}
Gene* search_for_i_gene(Genome *genome,int gene){
    struct Gene* pointer_gene = genome->head;
    for (int i=1; i<gene;i++){
        pointer_gene = pointer_gene->next;
    }
    return pointer_gene;
}
Gene* search_for_i_gene_previous(Genome *genome,int position){
    struct Gene* pointer_gene = genome->head;
    if(position>=2){
        for (int i=1; i<position-1;i++){
            pointer_gene = pointer_gene->next;
        }
    }
    return pointer_gene;
}
void delete_gene_in_genome(Genome *genome,int position){
    struct Gene *pointer_gene, *deleted_gene;
    if(position==1){
        deleted_gene = genome->head;
        genome->head = genome->head->next;
        free(deleted_gene);
    }else{
        pointer_gene = search_for_i_gene_previous(genome,position);
        deleted_gene = pointer_gene->next; 
        pointer_gene->next = pointer_gene->next->next;
        free(deleted_gene);
    }
    genome->size--;
}
void add_gene_to_genome(Genome *genome,Gene *gene, int position){
    struct Gene *pointer_gene, *next_connected;
    if(position==1){
        gene->next = genome->head;
        genome->head = gene; 
    }else{
        pointer_gene = search_for_i_gene_previous(genome,position);
        gene->next = pointer_gene->next;
        pointer_gene->next = gene;
    }
    genome->size++;
}
void mutation_change_weights(Genome * genome){
    struct Gene* gene = genome->head;
    float random_mutation; 
    for(int i=1; i<=genome->size; i++){
        random_mutation = ((float)rand()/(float)RAND_MAX);
        if(random_mutation <= chance_mutation_weight_per_genome){
            gene->weight += (((float)rand()/(float)RAND_MAX)*2*variation_weights)-variation_weights;
        }
        if(gene->is_bias){
            random_mutation = ((float)rand()/(float)RAND_MAX);
            if(random_mutation <= chance_mutation_weight_per_genome){
                gene->bias += (((float)rand()/(float)RAND_MAX)*2*variation_bias)-variation_bias;
            }
        }
        if(gene->is_activation){
            random_mutation = ((float)rand()/(float)RAND_MAX);
            if(random_mutation <= chance_mutation_weight_per_genome){
                gene->activation += (((float)rand()/(float)RAND_MAX)*2*variation_activation)-variation_activation;
                if(gene->activation > 0.5){
                    gene->activation = 0.5;
                }else if(gene->activation < 0){
                    gene->activation = 0;
                }
            }
        }
        gene = gene->next;
    }
}
void mutation_add_node(Genome *genome){
    int size = genome->size;
    struct Gene* pointer_gene;
    int gene =ceil(((float)rand()/(float)RAND_MAX)*size);
    pointer_gene = search_for_i_gene(genome,gene);
    int new_neuron = innovation+number_inputs+number_outputs;
    struct Gene* new_gene_1 = create_new_gene(pointer_gene->into,new_neuron,2,0,3,0,false,3,0,false,innovation,NULL);
    innovation++;
    struct Gene* new_gene_2 = create_new_gene(new_neuron, pointer_gene->out,1,pointer_gene->weight,1,pointer_gene->bias,pointer_gene->is_bias,1,pointer_gene->activation,pointer_gene->is_activation,innovation,NULL);
    innovation++;
    hitch_neuron_to_genome(new_neuron,genome);
    delete_gene_in_genome(genome,gene);
    add_gene_to_genome(genome,new_gene_2,gene);
    add_gene_to_genome(genome,new_gene_1,gene);
    genome->hiden_neurons++;
}
void get_random_final_output_position(Genome *genome, int excluded_neuron,Gene **left_gene, int *left_index,int *right_index, int *new_neuron){
    struct Gene *head_gene = genome->head;
    int size = genome->size;
    struct Gene *current_gene = head_gene;
    int left_in, right_in;
    struct Gene *left_ge;
    struct Neuron *neuron = genome->neurons;
    int i;
    int random_gene = (((float)rand()/(float)RAND_MAX)*genome->hiden_neurons)+1;
    if(random_gene == genome->hiden_neurons+1){
        random_gene--;
    }
    i=1;
    while(i<random_gene){
        neuron = neuron->next;
        i++;
    }
    if(neuron->identifier == excluded_neuron){
        if(neuron->next != NULL){
            neuron = neuron->next;
        }else{
            random_gene--;
            neuron = genome->neurons;
            i=1;
            while(i<random_gene){
                neuron = neuron->next;
                i++;
            }
        }
    }
    bool finded_right = false;
    i=1;
    while(!finded_right){
        if(current_gene->out == neuron->identifier){
            left_ge = current_gene;
            left_in = i;
        }
        if(current_gene->into == neuron->identifier){
            right_in = i;
            finded_right = true;
        }
        current_gene = current_gene->next;
        i++;
    }
    *left_gene = left_ge;
    *left_index = left_in;
    *right_index = right_in;
    *new_neuron = neuron->identifier;
}
bool check_new_conection(Gene *pointer_gene, int in, int out,int size, Gene *head_gene){
    bool finded = false;
    int i=1;
    struct Gene *current_gene = NULL;
    if(size >= 1){
        if(pointer_gene!=NULL){
            current_gene = pointer_gene->next;
        }else{
            current_gene = head_gene;
        }
    }
    while(!finded && i<=size){
        if(current_gene->into == in && current_gene->out == out){
            finded = true;
        }
        current_gene = current_gene->next;
        i++;
    }
    return finded;
}
void add_gene_to_genome_relative(Gene *pointer_gene,Gene *new_gene,int size, Genome *genome){
    int random_position;
    struct Gene *current_gene, *previous_gene;
    if(size >= 2){
        random_position = ((((float)rand()/(float)RAND_MAX))*size)+1;
    }else{
        random_position = size;
    }
    if(pointer_gene != NULL){
        current_gene = pointer_gene;
    }else{
        current_gene = genome->head;
    }
    if(random_position >=2){
        for(int i=1; i<random_position;i++){
            previous_gene = current_gene;
            current_gene = current_gene->next;
        }
    }else{
        if(pointer_gene != NULL){
            current_gene = pointer_gene->next;
            previous_gene = pointer_gene;
        }
    }
    if(pointer_gene != NULL){
        new_gene->next = current_gene;
        previous_gene->next = new_gene;
    }else{
        new_gene->next = genome->head;
        genome->head = new_gene;
    }
    genome->size++;
}
void mutation_add_connection(Genome *genome){  
    float random_type_connection = ((float)rand()/(float)RAND_MAX);
    int random_input_or_output = ((((float)rand()/(float)RAND_MAX))*2)+1;
    bool hiden_to_hiden_fail = false, hiden_to_external_fail = false, external_to_external_fail = false;
    struct Gene *gene_1 = NULL,*gene_2 = NULL, *left_gene = NULL;
    struct Gene *new_gene;
    int left_index_1,right_index_1, neuron_1; 
    int left_index_2, right_index_2, neuron_2;
    int relative_position, neuron_in, neuron_ou;

    if(random_type_connection <= structural_connections[0]){
        if(genome->hiden_neurons >= 2){
            get_random_final_output_position(genome,0,&gene_1,&left_index_1,&right_index_1,&neuron_1);
            get_random_final_output_position(genome,neuron_1,&gene_2,&left_index_2,&right_index_2,&neuron_2);
            if(left_index_1 < left_index_2){
                relative_position = right_index_2-left_index_1-1;
                neuron_in = neuron_1;
                neuron_ou = neuron_2;
                left_gene = gene_1;
            }else{
                relative_position = right_index_1-left_index_2-1;
                neuron_in = neuron_2;
                neuron_ou = neuron_1;
                left_gene = gene_2;
            }
            if(check_new_conection(left_gene,neuron_in,neuron_ou,relative_position,genome->head)){
                hiden_to_hiden_fail = true;
            }else{
                new_gene = create_new_gene(neuron_in,neuron_ou,3,0,2,0,false,2,0,false,innovation,NULL);
                innovation++;
                add_gene_to_genome_relative(left_gene,new_gene,relative_position,genome);
            }
        }else{
            hiden_to_hiden_fail = true;
        }
    }if((random_type_connection <= structural_connections[1] && random_type_connection > structural_connections[0]) || hiden_to_hiden_fail){
        if(genome->hiden_neurons >= 1){
            get_random_final_output_position(genome,0,&gene_1,&left_index_1,&right_index_1,&neuron_1);
            if(random_input_or_output == 1){
                relative_position = right_index_1-1;
                neuron_in = ((((float)rand()/(float)RAND_MAX))*number_inputs)+1;
                if(neuron_in==number_inputs+1){
                    neuron_in--;
                }
                neuron_ou = neuron_1;
                left_gene = NULL;
            }else if(random_input_or_output == 2){
                relative_position = genome->size-left_index_1;
                neuron_in = neuron_1;
                neuron_ou = ((((float)rand()/(float)RAND_MAX))*number_outputs)+1+number_inputs;
                if(neuron_ou==number_inputs+number_outputs+1){
                    neuron_ou--;
                }
                left_gene = gene_1;         
            }
            if(check_new_conection(left_gene,neuron_in,neuron_ou,relative_position,genome->head)){
                hiden_to_external_fail = true;
            }else{
                if(neuron_in>number_inputs && !genome->occupied_outputs[neuron_ou-number_inputs-1]){
                    genome->occupied_outputs[neuron_ou-number_inputs-1]=true;
                    new_gene = create_new_gene(neuron_in,neuron_ou,3,0,3,0,false,3,0,false,innovation,NULL);
                }else{ 
                    new_gene = create_new_gene(neuron_in,neuron_ou,3,0,2,0,false,2,0,false,innovation,NULL);
                }
                innovation++;
                add_gene_to_genome_relative(left_gene,new_gene,relative_position,genome);
            }
        }else{
            hiden_to_external_fail = true;
        }
    }if(random_type_connection > structural_connections[1] || hiden_to_external_fail){
        relative_position = genome->size;
        neuron_in = ((((float)rand()/(float)RAND_MAX))*number_inputs)+1;
        neuron_ou = ((((float)rand()/(float)RAND_MAX))*number_outputs)+1+number_inputs;
        if(neuron_ou==number_inputs+number_outputs+1){
            neuron_ou--;
        }
        if(neuron_in==number_inputs+1){
            neuron_in--;
        }
        left_gene = NULL;
        if(check_new_conection(left_gene,neuron_in,neuron_ou,relative_position,genome->head)){
            external_to_external_fail = true;
        }else{
            if(!genome->occupied_outputs[neuron_ou-number_inputs-1]){
                genome->occupied_outputs[neuron_ou-number_inputs-1]=true;
                new_gene = create_new_gene(neuron_in,neuron_ou,3,0,3,0,false,3,0,false,innovation,NULL);
            }else{
                new_gene = create_new_gene(neuron_in,neuron_ou,3,0,2,0,false,2,0,false,innovation,NULL);
            }
            innovation++;
            add_gene_to_genome_relative(left_gene,new_gene,relative_position,genome);
        }
    }
}
void delete_genome(Genome *genome){
    struct Gene *gene = genome->head, *delete_gene;
    // printf("heres 1\n");
    for(int i=1; i<=genome->size; i++){
        delete_gene = gene;
        gene = gene->next;
        free(delete_gene);
    }
    // printf("heres 2\n");
    delete_all_neurons(genome);
    // printf("heres 3\n");
    free(genome);
    kills++;
}
//-----------------------------------Funciones para la evaluacion de la red neuronal----------------------------
//-----------------------------------Funciones del arbol--------------------------------------------------------
Eval_tree* create_new_node_tree(int neuron,Eval_tree *parent_tree){
    struct Eval_tree *node = malloc(sizeof(Eval_tree));
    node->neuron = neuron;
    node->value = 0;
    node->activation = 0;
    node->parent_tree = parent_tree;
    node->left_tree = NULL;
    node->right_tree = NULL;
    return node;
}
Eval_tree* add_new_node_to_tree(Eval_tree* tree,int neuron){
    struct Eval_tree *pointer_tree = tree, *parent_tree;
    int right = 0;
    bool finded = false;
    if(tree == NULL){
        tree = create_new_node_tree(neuron,NULL);
    }else{
        while(!finded){
            if (pointer_tree == NULL){
                pointer_tree = create_new_node_tree(neuron,parent_tree);
                if(right==1){
                    parent_tree->right_tree = pointer_tree;
                }else{
                    parent_tree->left_tree = pointer_tree;
                }
                finded = true;
            }else{
                parent_tree = pointer_tree;
                if(pointer_tree->neuron > neuron){
                    pointer_tree = pointer_tree->left_tree;
                    right = 0;
                }else if(pointer_tree->neuron < neuron){
                    pointer_tree = pointer_tree->right_tree;
                    right = 1;
                }else{
                    finded = true;
                }
            }
        }
    }
    return tree;
}
void print_tree(Eval_tree *tree){
    printf(" Node:%d ",tree->neuron);
    printf(" (Left ->");
    if(tree->left_tree == NULL){
        printf(" NULL");
    }else{
        print_tree(tree->left_tree);
    }
    printf(")");
    printf("(Right ->");
    if(tree->right_tree == NULL){
        printf(" NULL");
    }else{
        print_tree(tree->right_tree);
    }
    printf(")");
}
Eval_tree* create_eval_tree(Genome *genome){
    struct Eval_tree *tree=NULL;
    int size = genome->size;
    int half = (size/2);
    Gene *pointer_gene = genome->head;
    for(int i=1; i<=half; i++){
        pointer_gene = pointer_gene->next;
    }
    for(int i=half; i<size;i++){
        tree = add_new_node_to_tree(tree, pointer_gene->into);
        tree = add_new_node_to_tree(tree, pointer_gene->out);
        pointer_gene = pointer_gene->next; 
    }
    pointer_gene = genome->head;
    for(int i=1; i<=half; i++){
        tree = add_new_node_to_tree(tree, pointer_gene->into);
        tree = add_new_node_to_tree(tree, pointer_gene->out);
        pointer_gene = pointer_gene->next; 
    }
    return tree;
}
void delete_tree(Eval_tree *tree){
    struct Eval_tree *pointer_tree=tree, *parent_tree=NULL;
    bool deleted=false;
    while(!deleted){
        if((pointer_tree->left_tree == NULL) && (pointer_tree->right_tree == NULL)){
            if(parent_tree==NULL){
                deleted=true;
            }else{
                if(pointer_tree->parent_tree->left_tree ==NULL){
                    pointer_tree->parent_tree->right_tree = NULL;
                }else{
                    pointer_tree->parent_tree->left_tree = NULL;
                }
            }
            free(pointer_tree);
            if(parent_tree!=NULL){
                pointer_tree = parent_tree;
                parent_tree = pointer_tree->parent_tree;
            }
        }else{
            parent_tree = pointer_tree;
            if(pointer_tree->left_tree != NULL){
                pointer_tree = pointer_tree->left_tree;
            }else if(pointer_tree->right_tree != NULL){
                pointer_tree = pointer_tree->right_tree;
            }
        }
    }
}
//--------------------------------------Evaluar red neuronal------------------------------------------------------
double see_value_in_tree(Eval_tree *tree,int neuron){
    double result=0;
    bool finded = false;
    Eval_tree *pointer_tree = tree;
    while((!finded) && pointer_tree!=NULL){
        if(pointer_tree->neuron == neuron){
            result = pointer_tree->value;
            finded = true;
        }else if(neuron < pointer_tree->neuron){
            pointer_tree = pointer_tree->left_tree;
        }else if(neuron > pointer_tree->neuron){
            pointer_tree = pointer_tree->right_tree;
        }
    }
    return result;
}
float see_activation_in_tree(Eval_tree *tree,int neuron){
    float result=0;
    bool finded = false;
    Eval_tree *pointer_tree = tree;
    while((!finded) && pointer_tree!=NULL){
        if(pointer_tree->neuron == neuron){
            result = pointer_tree->activation;
            finded = true;
        }else if(neuron < pointer_tree->neuron){
            pointer_tree = pointer_tree->left_tree;
        }else if(neuron > pointer_tree->neuron){
            pointer_tree = pointer_tree->right_tree;
        }
    }
    return result;
}
double add_value_in_tree(Eval_tree *tree,int neuron,double value,double activation){
    bool finded = false;
    double last_value;
    Eval_tree *pointer_tree = tree;
    while((!finded) && pointer_tree!=NULL){
        if(pointer_tree->neuron == neuron){
            last_value = pointer_tree->value;
            pointer_tree->value += value;
            pointer_tree->activation = activation;
            finded = true;
        }else if(neuron < pointer_tree->neuron){
            pointer_tree = pointer_tree->left_tree;
        }else if(neuron > pointer_tree->neuron){
            pointer_tree = pointer_tree->right_tree;
        }
    }
    return last_value;
}
double activation_function(double synapsis,float activation){
    float new_activation = round(activation*10)/10;
    float result;
    if(new_activation!=0){
        if(synapsis<(-1*new_activation)){
            result = 0;
        }else if(new_activation<synapsis){
            result = 1;
        }else{
            result = synapsis/(2*new_activation)+0.5;
        }
    }else{
        if(synapsis>0){
            result = 1;
        }else{
            result = 0;
        }
    }
    return result;
}
float* evaluate_network(Genome *genome, unsigned char (*inputs)[length][4],float *results,FILE *file_evaluation, int k){
    struct Eval_tree *tree = create_eval_tree(genome);
    struct Gene* gene = genome->head;
    double last_value;
    int coordinates[3];
    int neuron;
    // fprintf(file_evaluation,"The name is: %d and the image:%d\n",genome->name,k+1);
    for(int i=0; i<height; i++){
        for(int j=0; j<length; j++){
            for(int k=0; k<4; k++){
                neuron = i*(length*4) + j*4 + k+1;
                add_value_in_tree(tree,neuron,(double)inputs[i][j][k]/(double)255,0);
            }
        }
    }
    double synapsis;
    for(int i=1; i<=genome->size; i++){
        // fprintf(file_evaluation,"------------------------------------------------------------\n");
        // if(gene->into <= number_inputs){
        //     find_coordenates(gene->into,coordinates);
        //     fprintf(file_evaluation,"Gene:%d -> in:(%d,%d,%d)%d  out:%d\n",i,coordinates[0]+1,coordinates[1]+1,coordinates[2],gene->into,gene->out);
        // }else{
        //     fprintf(file_evaluation,"Gene:%d -> in:%d  out:%d\n",i,gene->into,gene->out);
        // }
        synapsis = see_value_in_tree(tree,gene->into);
        // fprintf(file_evaluation,"Input from neuron %d: %f",gene->into,synapsis);
        if(gene->into > number_inputs){
            synapsis = activation_function(synapsis,see_activation_in_tree(tree,gene->into));
            // fprintf(file_evaluation,"(act:%f) -> %f \n",see_activation_in_tree(tree,gene->into),synapsis);
        }else{
            // fprintf(file_evaluation,"-> %f \n",synapsis);
        }
        last_value = add_value_in_tree(tree,gene->out,synapsis*gene->weight,0);
        // fprintf(file_evaluation,"Neuron %d: %f + (%f * %f) = %f + %f = %f \n",gene->out,last_value,synapsis,gene->weight,last_value,synapsis*gene->weight,last_value+synapsis*gene->weight);
        if(gene->is_bias){
            last_value = add_value_in_tree(tree,gene->out,-1*gene->bias,0);
            // fprintf(file_evaluation,"bias: %f - %f = %f \n",last_value,gene->bias,last_value-gene->bias);
        }
        if(gene->is_activation){
            last_value = add_value_in_tree(tree,gene->out,0,gene->activation);
        }
        gene = gene->next;
    }
    // fprintf(file_evaluation,"------------------------------------------------------------\n");
    for(int i=0; i<number_outputs; i++){
        synapsis = see_value_in_tree(tree,number_inputs+i+1);
        // fprintf(file_evaluation,"output from neuron %d: %f (act:%f)-> ",number_inputs+i+1,synapsis,see_activation_in_tree(tree,number_inputs+i+1));
        synapsis = activation_function(synapsis,see_activation_in_tree(tree,number_inputs+i+1));
        // fprintf(file_evaluation,"%f \n",synapsis);
        *(results+i)= synapsis;
    }
    // fprintf(file_evaluation,"\n\n ");
    delete_tree(tree);
    return results;
}
//-----------------------------------Funciones fake para el fitness---------------------------------------
double fitness_increasing_function(int distance,int maximum_reward,int minimum_distance,int increasing){
    double fitness;
    if(distance > minimum_distance){
        fitness = 0;
    }else{
        fitness = ((double) maximum_reward / (double) pow(minimum_distance,increasing)) * pow((double)(distance-minimum_distance),increasing)*pow(-1,increasing);
    }
    return fitness;
}
double evaluate_false_fitness(unsigned char (*inputs)[length][4],int outputs, int intensities[2] ,Genome *genome,int (*sentence)[tests_numbers],int tes,FILE *file_evaluation, FILE *file_tests, int generation){
    float real_outputs[number_outputs];
    float result_angle;
    int choose_point, fitness = 0;
    int random_input = (((float)rand()/(float)RAND_MAX)*number_input_images);
    evaluate_network(genome,inputs,real_outputs,file_evaluation,random_input);
    result_angle = (float)((360*reaction_time)*(real_outputs[1]-real_outputs[0])*-max_velocity)/(float)(2*M_PI*car_length);
    if(abs(result_angle)<=max_angle){
        choose_point = (tan(result_angle*(M_PI/180))*(length/2)/tan(max_angle*(M_PI/180)))+length/2;
        if(choose_point==outputs){
            fitness = 100;
        }else{
            zeros++;
        }
    }
    sentence[0][tes]=outputs;
    sentence[1][tes]=choose_point;
    // if(generation>0){
    //     fprintf(file_tests,"Name: %d, input:%d %d,outputs: %.3f %.3f, velocity:%.2f %.2f, angle: %.1f, chose_points: %d, fitness:%d \n",genome->name,intensities[0],intensities[1], real_outputs[0],real_outputs[1],real_outputs[0]*max_velocity,real_outputs[1]*max_velocity,result_angle,choose_point,fitness);
    // }
    return fitness;
}
void create_false_fitness(unsigned char (*inputs)[length][4],int *output, int (*intensities),FILE *file_inputs, int name){
    int acu1 = 0, acu2 = 0;
    int rand1,rand2,res;
    int max_iter = 0;
    float color;
    // fprintf(file_inputs,"\n name: %d\n",name);
    for(int i=0; i<length; i++){
        intensities[i] = 0;
    }
    for(int i=0; i<height; i++){
        for(int j=0;j<length; j++){
            color = ((float)rand()/(float)RAND_MAX);
            if(color<color_propability){
                inputs[i][j][1] = ((float)rand()/(float)RAND_MAX)*255;                                  //red
                inputs[i][j][0] = ((float)rand()/(float)RAND_MAX)*inputs[i][j][1];                      //green
                inputs[i][j][0] = ((float)rand()/(float)RAND_MAX)*(inputs[i][j][1]-inputs[i][j][0]);    //blue
            }else{
                inputs[i][j][0] = ((float)rand()/(float)RAND_MAX)*255;          //red
                inputs[i][j][1] = ((float)rand()/(float)RAND_MAX)*255;          //green
                inputs[i][j][2] = ((float)rand()/(float)RAND_MAX)*255;          //blue
            }
            inputs[i][j][3] = j;
            
            // fprintf(file_inputs,"r:%d g:%d b:%d i:%d \t",inputs[i][j][0],inputs[i][j][1],inputs[i][j][2],inputs[i][j][3]);
            if(inputs[i][j][1] > (inputs[i][j][0]+inputs[i][j][2])){
                intensities[j] += inputs[i][j][1]-(inputs[i][j][0]+inputs[i][j][2]);
            }
        }
        // fprintf(file_inputs,"\n");
    }
    // fprintf(file_inputs,"intensities: %d %d\n",intensities[0],intensities[1]);
    for(int i=0; i<length; i++){
        if(intensities[i] > intensities[max_iter]){
            max_iter = i;
        }
    }
    *output = max_iter;
}
void check_sentence(int sentence_outputs[2][tests_numbers],Genome *genome){
    bool predic_out=true, real_out=true;
    int first_predic_out = sentence_outputs[1][0], first_real_out = sentence_outputs[0][0];
    int i=0;
    while(predic_out && i<tests_numbers){
        if(first_predic_out !=sentence_outputs[1][i]){
            predic_out=false;
        }
        i++;
    }
    i=0;
    while(real_out && i<tests_numbers){
        if(first_real_out !=sentence_outputs[0][i]){
            real_out=false;
        }
        i++;
    }
    if(predic_out && !real_out){
        genome->death_sentence=true;
        damned++;
    }
}
//-----------------------------------Funciones para el manejo de las especies----------------------------
//------------------------------------Eliminacion de individuos --------------------------------------------
void evaluate_fitness_of_all_species(Specie *specie, unsigned char (*inputs)[length][4], FILE *file_evaluation,FILE *file_tests,FILE *file_inputs,int generation){
    int outputs; 
    int sentence_outputs[2][tests_numbers];
    int intensities[length]; 
    int fitness, total_fitness;  
    struct Specie *current_specie = specie;
    damned = 0;
    while(current_specie!=NULL){
        struct Genome *current_genome = current_specie->first;
        int population = current_specie->population;
        for(int i=0; i<population; i++){
            total_fitness = 0;
            average_size += current_genome->size;
            for(int j=0; j<tests_numbers; j++){
                create_false_fitness(inputs,&outputs,intensities,file_inputs,current_genome->name);
                fitness = evaluate_false_fitness(inputs,outputs,intensities,current_genome,sentence_outputs,j,file_evaluation, file_tests,generation);
                total_fitness += fitness;
            }
            check_sentence(sentence_outputs,current_genome);
            current_genome->age++;
            current_genome->fitness = (double) total_fitness/(double)tests_numbers;
            current_genome = current_genome->next;
        }
        current_specie = current_specie->next;
    }
    // if(generation>0){
    //     fprintf(file_tests,"\n\n",fitness);
    // }   
}   
void search_position_order_genome(Genome **genome_head,Genome *genome){
    genome->next = NULL;
    genome->previous = NULL;
    bool finded = false;
    struct Genome *pointer_genome = *genome_head;
    struct Genome *previous_genome;
    if(pointer_genome == NULL){
        *genome_head = genome;
    }else{
        while ((pointer_genome!=NULL)&&(!finded)){
            if(pointer_genome->fitness <= genome->fitness){
                finded = true;
                genome->next = pointer_genome;
                genome->previous = pointer_genome->previous;
                pointer_genome->previous = genome;
                if(genome->previous!=NULL){
                    genome->previous->next = genome;
                }else{
                    *genome_head = genome;
                }
            }
            previous_genome = pointer_genome;
            pointer_genome = pointer_genome->next;
        }
        if(!finded){
            previous_genome->next = genome;
            genome->previous = previous_genome;
        }
    }
}
void order_by_fitness(Specie *specie){
    struct Specie *current_specie = specie;
    int killer_number, survivers;
    while(current_specie!=NULL){
        struct Genome *next_genome;
        struct Genome *current_genome = current_specie->first;
        struct Genome *previous_connected = current_specie->first->previous;
        struct Genome *next_connected = current_specie->last->next;
        struct Genome *new_first = NULL;
        int population = current_specie->population;
        for(int i=1; i<=population; i++){
            if(current_genome->previous!=NULL){
                current_genome->previous->next = current_genome->next;
            }
            if(current_genome->next!=NULL){
                current_genome->next->previous = current_genome->previous;
            }
            next_genome = current_genome->next;
            if(current_genome->death_sentence){
                delete_genome(current_genome);
                current_specie->population--;
            }else{
                search_position_order_genome(&new_first,current_genome);
            }
            current_genome = next_genome;
        }
        struct Genome *new_last = new_first;
        killer_number = floor((float)population*population_mortality);
        if(kills == population){
            printf(" toda la poblacion ha muerto en order by fitness\n");
        }
        if(damned > killer_number){
            survivers = population - damned;
        }else{
            survivers = population - killer_number;
        }
        population = current_specie->population;    
        for(int i=1; i<population; i++){
            new_last->rank = survivers-i+1;
            new_last = new_last->next;
        }
        current_specie->first = new_first;
        current_specie->last = new_last;
        new_first->previous = previous_connected;
        if(previous_connected!=NULL){
            previous_connected->next = new_first;
        }
        new_last->next = next_connected;
        if(next_connected!=NULL){
            next_connected->previous = new_last;
        }
        current_specie = current_specie->next;
    }
}
void kill_individuals(Specie *specie){
    struct Specie *current_specie = specie;
    while(current_specie != NULL){
        struct Genome *current_genome = current_specie->first,*next_genome;
        struct Genome *next_specie_genome = current_specie->last->next;
        int population = current_specie->population;
        int killer_number;
        int killed_percentage = floor((float)(population+damned)*population_mortality);
        int living_number = population+damned-killed_percentage-1;
        current_specie->reproductivity = 0;
        if(damned >= killed_percentage){
            living_number = population-1;
        }else{
            killer_number = killed_percentage-kills;
        }
        for(int i=0; i<living_number;i++){
            current_specie->reproductivity += current_genome->rank + (age_importance*current_genome->age); 
            current_genome = current_genome->next;
        }
        current_specie->reproductivity += current_genome->rank + (age_importance*current_genome->age);
        current_specie->last = current_genome;
        if(damned < killed_percentage){
            current_genome = current_genome->next;
            for(int i=0; i<killer_number;i++){
                next_genome = current_genome->next;
                delete_genome(current_genome);
                current_genome = next_genome;
            }
            current_specie->population = population-killer_number;
        }        
        current_specie->last->next = next_specie_genome;
        if(next_specie_genome!=NULL){
            next_specie_genome->previous = current_specie->last;
        }
        current_specie = current_specie->next;
    }
}
//------------------------------------Reproduccion de especies ---------------------------------------------
void calculate_fitness_of_species(Specie *specie){
    struct Specie *current_specie = specie;
    total_fitness = 0;
    while(current_specie != NULL){
        struct Genome *current_genome = current_specie->first;
        int population = current_specie->population;
        current_specie->fitness = 0;
        for(int i=0; i<population; i++){
            current_genome->sons = 0;
            current_specie->fitness += current_genome->fitness/population;
            total_fitness += current_genome->fitness/population;
            current_genome = current_genome->next;
        }
        current_specie = current_specie->next;
    }
}
void calculate_offspring(Specie *specie){
    struct Specie *current_specie = specie;   
    double ratio = total_fitness/(double)(kills);
    if(ratio == 0){
        ratio = (double)kills/(double)reproductive_species;
    }
    int new_offspring = 0;
    while(current_specie != NULL){
        if(total_fitness!=0){
            current_specie->reproductions = floor((current_specie->fitness)/(ratio/2));
        }else{
            current_specie->reproductions = ratio;
        }
        new_offspring += current_specie->reproductions;
        current_specie->new_childs = 0;
        current_specie = current_specie->next;
    }
    current_specie = specie;
    for(int i=1; i<=(kills-new_offspring); i++){
        current_specie->reproductions++;
        current_specie = current_specie->next;
    }
}
//---------------------------------------Selecting partner---------------------------------------------------
void choose_genome(Specie *specie, Genome **partner1, Genome **partner2){
    struct Genome *current_genome = specie->first;
    int random_reproductivity_1, random_reproductivity_2;
    int iter_reproductivity = 0;
    bool finded_1 = false, finded_2 = false;
    random_reproductivity_1 = (((float)rand()/(float)RAND_MAX)*specie->reproductivity)+1;
    random_reproductivity_2 = (((float)rand()/(float)RAND_MAX)*specie->reproductivity)+1;

    while(!finded_1 || !finded_2){
        iter_reproductivity += current_genome->rank +(age_importance*current_genome->age);
        if(random_reproductivity_1 <= iter_reproductivity && !finded_1){
            *partner1 = current_genome;
            finded_1 = true;
        }else if(random_reproductivity_2 <= iter_reproductivity && !finded_2){
            *partner2 = current_genome;
            finded_2 = true;
        }else{
            current_genome = current_genome->next;
        }
    }
}
void select_partner(Specie *current_specie, Genome **parent1, Genome **parent2){
    current_specie->reproductions--;
    current_specie->reproductions--;

    choose_genome(current_specie,parent1,parent2);
    (*parent1)->sons++;
    (*parent2)->sons++;
    // printf("parent: s = %d n = %d,            ",current_specie->identifier,current_genome->name);
}
//--------------------------------------Sexual reproduction--------------------------------------------------
bool is_in_tree(int neuron,Eval_tree *tree){
    struct Eval_tree *pointer_tree = tree;
    bool finded = false;
    while(pointer_tree!=NULL && (!finded)){
        if(pointer_tree->neuron == neuron){
            finded = true;
        }else{
            if(neuron < pointer_tree->neuron){
                pointer_tree = pointer_tree->left_tree;
            }else if(pointer_tree->neuron < neuron){
                pointer_tree = pointer_tree->right_tree;
            }
        }
    }
    return finded;
}
Gene* copy_gene(Gene *gene){
    struct Gene *new_gene = malloc(sizeof(Gene));
    new_gene->into = gene->into;
    new_gene->out = gene->out;
    new_gene->origin = gene->origin;
    new_gene->weight = gene->weight;
    new_gene->next = NULL;
    new_gene->bias = gene->bias;
    new_gene->is_bias = gene->is_bias;
    new_gene->activation = gene->activation;
    new_gene->is_activation = gene->is_activation;
    return new_gene;
}
Genome* copy_genome(Genome *parent){
    struct Genome *new_child = create_new_genome();
    struct Gene *copied_gene, *previous_copied_gene, *current_gene;
    current_gene = parent->head;
    new_child->size = parent->size;
    new_child->hiden_neurons = parent->hiden_neurons;
    new_child->occupied_outputs[0] = parent->occupied_outputs[0];
    new_child->occupied_outputs[1] = parent->occupied_outputs[1];
    for(int i=1; i<=parent->size; i++){
        copied_gene = copy_gene(current_gene);
        if(i==1){
            new_child->head = copied_gene;
        }else{
            previous_copied_gene->next = copied_gene;
        }
        previous_copied_gene = copied_gene;
        current_gene = current_gene->next;
    }
    copy_all_neurons(parent,new_child);
    return new_child;
}
void alter_genome(Genome *child, Genome *parent_2){
    struct Gene *current_gene = child->head;
    struct Gene *current_gene_parent;
    float probability;
    for(int i=1; i<=child->size; i++){
        current_gene_parent = parent_2->head;
        for(int j=1; j<=parent_2->size; j++){
            if(current_gene->origin == current_gene_parent->origin){
                probability = ((float)rand()/(float)RAND_MAX);
                if(probability < 0.5){
                    current_gene->weight = current_gene_parent->weight;
                    current_gene->is_bias = current_gene_parent->is_bias;
                    current_gene->bias = current_gene_parent->bias;
                    current_gene->is_activation = current_gene_parent->is_activation;
                    current_gene->activation = current_gene_parent->activation;
                }
            }
            current_gene_parent = current_gene_parent->next;
        }
        current_gene = current_gene->next;
    }
}
void chain_genes(Gene *gene1,int segment1, Gene *gene2,int segment2,Gene **head,Gene **tail){
    struct Gene *current_gene1 = gene1, *current_gene2 = gene2;
    struct Gene *copy_gene1, *copy_gene2;
    struct Gene *head1=NULL, *tail1=NULL, *head2=NULL, *tail2=NULL;
    *head = NULL; *tail = NULL;
    for(int i=1; i<segment1; i++){
        if(i==1){
            copy_gene1 = copy_gene(current_gene1);
            head1 = copy_gene1;
        }else{
            copy_gene1->next = copy_gene(current_gene1);
            copy_gene1 = copy_gene1->next;
        }
        if(i==segment1-1){
            tail1 = copy_gene1;
        }
        current_gene1 = current_gene1->next;
    }
    for(int i=1; i<segment2; i++){
        if(i==1){
            copy_gene2 = copy_gene(current_gene2);
            head2 = copy_gene2;
        }else{
            copy_gene2->next = copy_gene(current_gene2);
            copy_gene2 = copy_gene2->next;
        }
        if(i==segment2-1){
            tail2 = copy_gene2;
        }
        current_gene2 = current_gene2->next;
    }
    if(head1!=NULL && head2!=NULL){
        tail1->next = head2;
        *head = head1;
        *tail = tail2;
    }if(head1==NULL && head2!=NULL){
        *head = head2;
        *tail = tail2; 
    }
    if(head1!=NULL && head2==NULL){
        *head = head1;
        *tail = tail1; 
    }
}
void chain_final_genes(Gene *gene1,Gene *gene2,Gene **head,Gene **tail){
    struct Gene *current_gene1 = gene1, *current_gene2 = gene2;
    struct Gene *copy_gene1, *copy_gene2;
    struct Gene *head1=NULL, *tail1=NULL, *head2=NULL, *tail2=NULL;
    *head = NULL; *tail = NULL;
    
    while(current_gene1!=NULL){
        if(head1==NULL){
            copy_gene1 = copy_gene(current_gene1);
            head1 = copy_gene1;
        }else{
            copy_gene1->next = copy_gene(current_gene1);
            copy_gene1 = copy_gene1->next; 
        }
        tail1 = copy_gene1;
        current_gene1 = current_gene1->next;
    }
    while(current_gene2!=NULL){
        if(head2==NULL){
            copy_gene2 = copy_gene(current_gene2);
            head2 = copy_gene2;
        }else{
            copy_gene2->next = copy_gene(current_gene2);
            copy_gene2 = copy_gene2->next; 
        }
        tail2 = copy_gene2;
        current_gene2 = current_gene2->next;
    }
    if(head1!=NULL && head2!=NULL){
        tail1->next = head2;
        *head = head1;
        *tail = tail2;
    }if(head1==NULL && head2!=NULL){
        *head = head2;
        *tail = tail2; 
    }
    if(head1!=NULL && head2==NULL){
        *head = head1;
        *tail = tail1; 
    }
}
Genome* mix_genomes(Genome *genome1, Genome *genome2){
    struct Genome *new_child = create_new_genome();
    struct Eval_tree *tree_2 = create_eval_tree(genome2);
    int neuron_1;
    struct Gene *search_point1_genome1=NULL, *search_point2_genome1;
    struct Gene *search_point1_genome2=NULL, *search_point2_genome2;
    struct Gene *head_child, *body_child;
    struct Gene *head, *tail;
    float chance_add_node, chance_add_connection, chance_weight;
    bool search, search_for_match = true;
    int segment1 = 1, segment2 = 1,join=0;
    search_point2_genome1 = genome1->head;

    while(search_for_match){
        neuron_1 = search_point2_genome1->out;
        if(is_in_tree(neuron_1,tree_2)){
            if(search_point1_genome2 == NULL){
                search_point2_genome2 = genome2->head;
            }else{
                search_point2_genome2 = search_point1_genome2->next;
                segment2 = 1;
            }
            search = true;
            while(search){
                if(search_point2_genome1->origin==search_point2_genome2->origin){
                    if((search_point1_genome1==NULL)&&(search_point1_genome2==NULL)){
                        chain_genes(genome1->head,segment1, genome2->head,segment2,&head,&tail);
                        if((head != NULL) || (tail != NULL)){
                            head_child = head;
                            body_child = tail;
                            if(genome1->fitness > genome2->fitness){
                                body_child->next = copy_gene(search_point2_genome1);
                            }else{
                                body_child->next = copy_gene(search_point2_genome2);
                            }
                            body_child = body_child->next;
                        }else{
                            if(genome1->fitness > genome2->fitness){
                                head_child = copy_gene(search_point2_genome1);
                            }else{
                                head_child = copy_gene(search_point2_genome2);
                            }
                            body_child = head_child;
                        }
                    }else{
                        chain_genes(search_point1_genome1->next,segment1, search_point1_genome2->next,segment2,&head,&tail);
                        if((head != NULL) || (tail != NULL)){
                            body_child->next = head;
                            body_child = tail;
                            if(genome1->fitness > genome2->fitness){
                                body_child->next = copy_gene(search_point2_genome1);
                            }else{
                                body_child->next = copy_gene(search_point2_genome2);
                            }
                            body_child = body_child->next;
                        }else{
                            if(genome1->fitness > genome2->fitness){
                                body_child->next = copy_gene(search_point2_genome1);
                            }else{
                                body_child->next = copy_gene(search_point2_genome2);
                            }
                            body_child = body_child->next;
                        }
                    }
                    search = false;
                    search_point1_genome1 = search_point2_genome1;
                    search_point1_genome2 = search_point2_genome2;
                    segment1 = 0;
                    segment2 = 0;
                    join++;
                }else if(search_point2_genome2->next == NULL || search_point2_genome2->into == neuron_1){
                    search = false;
                    segment2 = 1;
                }else{
                    search_point2_genome2 = search_point2_genome2->next;
                    segment2++;
                }
            }

        }
        search_point2_genome1 = search_point2_genome1->next;
        segment1++;
        if(search_point1_genome2!=NULL){
            if(search_point1_genome2->next==NULL || search_point2_genome1==NULL){
                search_for_match = false;
            }
        }else{
            if(search_point2_genome1==NULL){
                search_for_match = false;
            }
        }
    }
    if(search_point1_genome1==NULL && search_point1_genome2==NULL){
        chain_final_genes(genome1->head,genome2->head,&head,&tail);
        head_child = head;
        body_child = tail;
    }else{
        chain_final_genes(search_point1_genome1->next,search_point1_genome2->next,&head,&tail);
        if((head != NULL) || (tail != NULL)){
            body_child->next = head;
            body_child = tail;
        }
    }
    delete_tree(tree_2);

    new_child->head = head_child;
    new_child->size = genome1->size + genome2->size - join;
    new_child->hiden_neurons = genome1->hiden_neurons + genome2->hiden_neurons;
    if(new_child->size > N){
        N = new_child->size;
    }
    copy_all_neurons_for_mix_genomes(genome1,genome2,new_child);
    return new_child;
}
Genome* sexual_reproduction(Genome *parent_1, Genome *parent_2){
    struct Genome *new_child;
    float chance_add_node, chance_add_connection, chance_weight;
    if (abs(parent_1->fitness-parent_2->fitness) > epsilon){
        if(parent_1->fitness > parent_2->fitness){
            new_child = copy_genome(parent_1);
            alter_genome(new_child,parent_2);
        }else if(parent_2->fitness > parent_1->fitness){
            new_child = copy_genome(parent_2);
            alter_genome(new_child,parent_1);
        }
        // if(new_child->name==51){
        //     printf("p1:%d + p2:%d -> c:%d\n",parent_1->name,parent_2->name,new_child->name);
        //     printf("no mixing\n");
        // }
    }else{
        new_child = mix_genomes(parent_1,parent_2);
        // if(new_child->name==51){
        //     printf("p1:%d + p2:%d -> c:%d\n",parent_1->name,parent_2->name,new_child->name);
        //     printf("mixing\n");
        // }
    }
    chance_add_node = ((float)rand()/(float)RAND_MAX);
    chance_add_connection = ((float)rand()/(float)RAND_MAX);
    chance_weight = ((float)rand()/(float)RAND_MAX);
    if(chance_add_node < chance_mutation_add_node){
        mutation_add_node(new_child);
    }
    if(chance_add_connection < chance_mutation_add_connection){
        mutation_add_connection(new_child);
    }
    if(chance_weight < chance_mutation_change_weights){
        mutation_change_weights(new_child);
    }
    return new_child;
}
//--------------------------------------Add child to population--------------------------------------------
float calculate_distance_between_genomes(Genome *genome1, Genome *genome2){
    struct Eval_tree *tree_2 = create_eval_tree(genome2);
    int neuron_1;
    struct Gene *search_point1_genome1 = NULL, *search_point2_genome1;
    struct Gene *search_point1_genome2 = NULL, *search_point2_genome2;
    int exces_genome1, exces_genome2, join=0;
    float genome_weight = 0.0;
    bool search, search_for_match = true;
    float distance;
    search_point2_genome1 = genome1->head;
    while(search_for_match){
        neuron_1 = search_point2_genome1->out;
        if(is_in_tree(neuron_1,tree_2)){
           if(search_point1_genome2 == NULL){
                search_point2_genome2 = genome2->head;
            }else{
                search_point2_genome2 = search_point1_genome2->next;
            }
            search = true;
            while(search){
                if(search_point2_genome1->origin==search_point2_genome2->origin){
                    search = false;
                    search_point1_genome1 = search_point2_genome1;
                    search_point1_genome2 = search_point2_genome2;
                    genome_weight = fabs(search_point2_genome1->weight - search_point2_genome2->weight);
                    join++;
                }else if(search_point2_genome2->next == NULL || search_point2_genome2->into == neuron_1){
                    search = false;
                }else{
                    search_point2_genome2 = search_point2_genome2->next;
                }
            }
        }
        search_point2_genome1 = search_point2_genome1->next;
        if(search_point1_genome2!=NULL){
            if(search_point1_genome2->next==NULL || search_point2_genome1==NULL){
                search_for_match = false;
            }
        }else{
            if(search_point2_genome1==NULL){
                search_for_match = false;
            }
        }
    }
    delete_tree(tree_2);
    exces_genome1 = genome1->size - join;
    exces_genome2 = genome2->size - join;
    if(join==0){
        distance = c1*(exces_genome1+exces_genome2);
        // distance = c1*((float)(exces_genome1+exces_genome2)/(float)N);
    }else{
        distance = (c1*(exces_genome1+exces_genome2)+c2*(genome_weight/join));
        // distance = (c1*((float)(exces_genome1+exces_genome2)/(float)N)+c2*(genome_weight/join));
    }
    return distance;
}
void add_distance_to_species(float distance, int id, Specie *specie){
    struct Distances *current_distance = specie->distances, *previous_distance=NULL;
    struct Distances *node_distance;
    bool finded = false;
    while((!finded) && (current_distance!=NULL)){
        if(distance<current_distance->distance){
            struct Distances *node_distance = malloc(sizeof(Distances));
            node_distance->identifier = id;
            node_distance->distance = distance;
            node_distance->next = current_distance;
            if(previous_distance!=NULL){
                previous_distance->next = node_distance;
            }else{
                specie->distances = node_distance;
            }
            finded = true;
        }
        previous_distance = current_distance;
        current_distance = current_distance->next;
    }if(!finded){
        struct Distances *node_distance = malloc(sizeof(Distances));
        node_distance->identifier = id;
        node_distance->distance = distance;
        node_distance->next = NULL;
        if(previous_distance!=NULL){
            previous_distance->next = node_distance;
        }else{
            specie->distances = node_distance;
        }
    }
}
void add_child_to_specie(Genome *genome, Specie *current_specie, Specie *previous_specie){
    genome->next = current_specie->first;
    if(previous_specie!=NULL){
        previous_specie->last->next = genome;
        genome->previous = previous_specie->last;
    }
    current_specie->first->previous = genome;
    current_specie->first = genome;
    // genome->name = current_specie->max_name;
    // current_specie->max_name++;
    maxim_name++;

    current_specie->new_childs++;
    current_specie->population++;
}
void create_new_specie(Genome *genome, Specie **head_specie){
    struct Specie *new_specie = malloc(sizeof(Specie));
    struct Specie *current_specie = *head_specie;
    struct Specie *previous_specie = NULL;
    float distance;
    new_specie->identifier = number_species;
    new_specie->population = 1;
    new_specie->fitness = 0;
    new_specie->reproductions = 0;
    new_specie->new_childs = 0;
    new_specie->first = genome;
    new_specie->last = genome;
    new_specie->previous_max_fitness = 0;
    new_specie->max_fitness = 0;
    new_specie->it_advance = 0;
    new_specie->max_name = 1;
    new_specie->next = NULL;
    new_specie->distances = NULL;
    add_distance_to_species(0,new_specie->identifier,new_specie);
    while(current_specie!=NULL){
        distance = calculate_distance_between_genomes(genome,current_specie->first);
        add_distance_to_species(distance,current_specie->identifier,new_specie);
        add_distance_to_species(distance,new_specie->identifier,current_specie);
        previous_specie = current_specie;
        current_specie = current_specie->next;
    }
    if(previous_specie!=NULL){
        previous_specie->next = new_specie;
        previous_specie->last->next = genome;
        genome->previous = previous_specie->last;
    }else{
        *head_specie = new_specie;
    }
    genome->name = maxim_name;
    maxim_name++;
    // genome->name = new_specie->max_name;
    // new_specie->max_name++;

    // printf("child: s = %d n = %d \n",new_specie->identifier,genome->name);
    number_species++;
    living_species++;
}
void add_child_to_population(Genome *genome, Specie **head_specie){
    struct Specie *current_specie = *head_specie, *previous_specie = NULL;
    float distance;
    bool finded=false;
    while((!finded) && (current_specie!=NULL)){
        distance = calculate_distance_between_genomes(genome,current_specie->first);
        if(distance < distance_compatibility_in_specie){
            add_child_to_specie(genome,current_specie,previous_specie);
            finded = true;
        }
        previous_specie = current_specie;
        current_specie = current_specie->next;
    }
    if(!finded){
        create_new_specie(genome,head_specie);
    }
}
//---------------------------------------Check stop species--------------------------------------------------
Specie* delete_specie(Specie *prev_specie,Specie *specie){
    struct Genome *current_genome = specie->first, *deleted_genome;
    struct Genome *first_genome = NULL, *last_genome = NULL;
    struct Specie *result_specie = NULL;
    struct Distances *dist, *deleted_distance;

    dist = specie->distances;
    while(dist!=NULL){
        deleted_distance = dist;
        dist = dist->next;
        free(deleted_distance);
    }
    //borrar distancias
    for(int i=1; i<=specie->population; i++){
        deleted_genome = current_genome;
        current_genome = current_genome->next;
        delete_genome(deleted_genome);
    }
    if(prev_specie != NULL){
        first_genome = prev_specie->last;
    }
    if(specie->next != NULL){
        last_genome = specie->next->first;
    }
    if(prev_specie != NULL){
        prev_specie->next = specie->next;
        first_genome->next = last_genome;
        if(last_genome != NULL){
            last_genome->previous = first_genome;
        }
    }else{
        result_specie = specie->next;
        if(last_genome != NULL){
            last_genome->previous = NULL;
        }
    }
    free(specie);
    living_species = living_species - 1;
    return result_specie;
}
void delete_distance(Specie *head_specie, int identifier){
    struct Specie *specie = head_specie;
    bool finded;
    struct Distances *dist, *deleted_dist, *prev_dist;
    while(specie!=NULL){
        if(specie->identifier != identifier){
            finded = false;
            dist = specie->distances;
            while(!finded){
                if(dist->identifier == identifier){
                    finded = true;
                    deleted_dist = dist;
                    prev_dist->next = dist->next;
                }
                prev_dist = dist;
                dist = dist->next;
            }
            free(deleted_dist);
        }
        specie = specie->next;
    }
}
void check_stop_species(Specie **head_specie){
    struct Specie *current_specie = *head_specie, *prev_specie = NULL;
    struct Specie *deleted_specie, *new_head;
    while(current_specie!=NULL){
        if(current_specie->it_advance == 0){
            current_specie->max_fitness = current_specie->fitness;
            current_specie->it_advance++;
            prev_specie = current_specie;
            current_specie = current_specie->next;
        }else{
            if(current_specie->fitness > current_specie->max_fitness){
                current_specie->max_fitness = current_specie->fitness;
            }
            if(current_specie->it_advance == generation_stop){
                if((current_specie->max_fitness - current_specie->previous_max_fitness) > advance_in_fitness){
                    current_specie->it_advance = 1;
                    current_specie->previous_max_fitness =  current_specie->max_fitness;
                    current_specie->max_fitness = current_specie->fitness;
                    prev_specie = current_specie;
                    current_specie = current_specie->next;
                }else{
                    deleted_specie = current_specie;
                    total_fitness = total_fitness - current_specie->fitness;
                    delete_distance(*head_specie, current_specie->identifier);
                    new_head = delete_specie(prev_specie,deleted_specie);
                    if(new_head!=NULL || current_specie->next == NULL){
                        *head_specie = new_head;
                    }
                    current_specie = current_specie->next;
                }
            }else{
                current_specie->it_advance++;
                prev_specie = current_specie;
                current_specie = current_specie->next;
            }
        }    
    }
}
//--------------------------------------Auxiliar functions----------------------------------------------------
void first_population(Specie **specie,int population,FILE *file_genomes){
    struct Genome *genome;
    for(int i=0; i<population; i++){
        genome = create_new_genome();
        mutation_add_connection(genome);
        mutation_add_connection(genome);
        mutation_add_connection(genome);
        add_child_to_population(genome,specie);
        // print_genome_in_files(genome, file_genomes);
    }
}
void print_species_in_files(Specie *specie,FILE *file_resume){
    struct Genome *first_genome;
    struct Genome *last_genome;
    struct Genome *current_genome;
    struct Distances *dist;
    while(specie!=NULL){
        // printf("----------------------------------------------------------------------------------\n");
        // printf("specie %d has population %d and fitness %f and reproductions %d and new childs %d and max fitness: %f, previous max: %f and iter: %d\n",specie->identifier, specie->population, specie->fitness,specie->reproductions,specie->new_childs,specie->max_fitness,specie->previous_max_fitness,specie->it_advance);
        dist = specie->distances;
        // printf("distances: ");
        // while(dist!=NULL){
        //     printf(" id = %d, dist = %f ->",dist->identifier,dist->distance);
        //     dist = dist->next;
        // }
        // printf("----------------------------------------------------------------------------------\n");
        first_genome = specie->first;
        last_genome = specie->last;
        current_genome = first_genome;
        for(int i=1; i<=specie->population; i++){
            if(i==1){
                // fprintf(file_resume,"The name is %d and the fitness is: %d and age: %d and reproductivity: %d \n",current_genome->name,(int)current_genome->fitness,current_genome->age,current_genome->rank + current_genome->age);
                fprintf(file_resume,"The name is %d and the fitness is: %d ",current_genome->name,(int)current_genome->fitness);            
            }
            // print_genome(current_genome);
            // print_genome_genes(current_genome,current_genome->size);
            current_genome = current_genome->next;
        }
        // printf("----------------------------------------------------------------------------------\n");
        specie = specie->next;
    }
}
//---------------------------------------Offspring functions----------------------------------------------------
void print_offspring_first(int (*offspring)[number_offspring]){
    for(int i=0;i<(1-population_mortality)*population;i++){
        offspring[i][0]=0;
    }
}
void print_offspring_second(int (*offspring)[number_offspring],int parent_name,int child_name){
    bool finded = false;
    int i=0;
    int brothers;
    while(!finded && offspring[i][0]!=0){
        if(offspring[i][0]==parent_name){
            brothers = offspring[i][1];
            offspring[i][brothers+2]=child_name;
            offspring[i][1]++;
            finded=true;
        }
        i++;
    }
    if(!finded){
        offspring[i][0]= parent_name;
        offspring[i][1]=1;
        offspring[i][2]= child_name;
    }
}
void print_offspring_third(int (*offspring)[number_offspring],FILE *file_offspring){
    int i=0;
    int brothers;
    while(offspring[i][0]!=0){
        brothers = offspring[i][1];
        fprintf(file_offspring,"%d: ",offspring[i][0]);
        for(int j=0;j<brothers;j++){
            fprintf(file_offspring,"%d,",offspring[i][j+2]);
        }
        fprintf(file_offspring,"\n");
        i++;
    }
}
int main(){
    srand(time(NULL));
    struct Specie *specie = NULL;
    struct Specie *current_specie;
    struct Genome *new_child;
    struct Genome *partner1, *partner2;
    int generation = 1;
    int offspring[50][number_offspring];


    // unsigned char (*inputs)[height][length][3] = malloc(number_input_images*height*length*3*sizeof(char));
    unsigned char (*inputs)[length][4] = malloc(height*length*4*sizeof(char));
    // int (*outputs)[length]= malloc(number_input_images*length*sizeof(int));

    FILE *file_genomes = fopen("genome.txt","w");                   //first_population(), bucle while   
    FILE *file_resume = fopen("resume.txt","w");                    //
    FILE *file_evaluation = fopen("evaluation.txt","w");            //evaluate_network()
    FILE *file_tests = fopen("tests.txt","w");                      //evaluate_fitness_of_all_species(), evaluate_false_fitness()
    FILE *file_offspring = fopen("offspring.txt","w");              //activar las 3 funciones print_offspring()
    FILE *file_inputs = fopen("inputs.txt","w");                    //bucle while, create_false_fitness()

    // charge_inputs_and_outputs(inputs,outputs);
    first_population(&specie,population,file_genomes);
    while(generation<=max_generations && (specie!=NULL)){
        start = clock();
        kills = 0;
        zeros = 0;
        average_size = 0;
        // fprintf(file_inputs,"\nGeneration %d: ",generation);
        evaluate_fitness_of_all_species(specie,inputs,file_evaluation,file_tests,file_inputs,generation);
        printf("Generacion: %d, zeros: %d/%d\n",generation,zeros,population*tests_numbers);
        order_by_fitness(specie);
        kill_individuals(specie);
        calculate_fitness_of_species(specie);
        check_stop_species(&specie);
        fprintf(file_resume,"Generation %d: ",generation);
        print_species_in_files(specie,file_resume);
        if(specie!=NULL){
            print_genome_in_files(specie->first, file_genomes);
            reproductive_species = living_species;
            calculate_offspring(specie);
            current_specie = specie;
            new_childs = 0;
            // fprintf(file_offspring,"Generation %d: \n",generation);
            // print_offspring_first(offspring);
            while(new_childs!=kills){
                if(current_specie->reproductions>0){
                    select_partner(current_specie, &partner1, &partner2);
                    new_child = sexual_reproduction(partner1,partner2);
                    add_child_to_population(new_child, &specie);
                    // print_genome_in_files(new_child, file_genomes);
                    // print_offspring_second(offspring,partner1->name,new_child->name);
                    new_childs++;
                }else{
                    current_specie = current_specie->next;
                }
            }
            // print_offspring_third(offspring,file_offspring);
            // fprintf(file_offspring,"\n\n");
            generation++;
            end = clock();
            cpu_time_used = ((double)(end-start))/CLOCKS_PER_SEC;
            fprintf(file_resume,"time: %f, size: %f \n",cpu_time_used,(double)average_size/population);
        }else{
            printf("todas las especies se han extinto\n");
        }
    }
    // free(inputs);
    // free(outputs);

    fclose(file_genomes);
    fclose(file_resume);
    fclose(file_inputs);
    fclose(file_offspring);
    fclose(file_evaluation);
    fclose(file_tests);
}

// gcc neat_asexual.c -o neat -lm
// gcc neat_sexual.c -o neat -lm