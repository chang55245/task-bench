/* Copyright 2020 Los Alamos National Laboratory
 * Licensed under the Apache License, Version 2.0 */

#include <stdarg.h>
#include <assert.h>
#include <string.h>
#include <algorithm> 
#include <unistd.h>
#include "core.h"
#include "timer.h"
// #include "dash.h"

#define VERBOSE_LEVEL 0
#define USE_CORE_VERIFICATION

extern "C" {
typedef struct tile_s {
  float dep;
  char *output_buff;
}tile_t;

typedef struct payload_s {
  int x;
  int y;
  TaskGraph graph;
}payload_t;

typedef struct matrix_s {
  tile_t *data;
  int M;
  int N;
}matrix_t;

// Single scratch memory buffer for serial version
char *scratch_memory = NULL;

static inline void task1(tile_t *tile_out, payload_t payload)
{
  if (tile_out == NULL || tile_out->output_buff == NULL) {
    fprintf(stderr, "Error: NULL pointer in task1 at position (%d, %d)\n", payload.x, payload.y);
    return;
  }

#if defined (USE_CORE_VERIFICATION)    
  TaskGraph graph = payload.graph;
  char *output_ptr = (char*)tile_out->output_buff;
  size_t output_bytes = graph.output_bytes_per_task;
  std::vector<const char *> input_ptrs;
  std::vector<size_t> input_bytes;
  input_ptrs.push_back((char*)tile_out->output_buff);
  input_bytes.push_back(graph.output_bytes_per_task);
  
  graph.execute_point(payload.y, payload.x, output_ptr, output_bytes,
                     input_ptrs.data(), input_bytes.data(), input_ptrs.size(), 
                     scratch_memory, graph.scratch_bytes_per_task);
#else  
  tile_out->dep = 0;
  printf("Task1 x %d, y %d, out %f\n", payload.x, payload.y, tile_out->dep);
#endif  
}

static inline void task2(tile_t *tile_out, const std::vector<tile_t*> &tile_in, payload_t payload)
{
  if (tile_out == NULL || tile_out->output_buff == NULL) {
    fprintf(stderr, "Error: NULL pointer in task2 output at position (%d, %d)\n", payload.x, payload.y);
    return;
  }

#if defined (USE_CORE_VERIFICATION)    
  TaskGraph graph = payload.graph;
  char *output_ptr = (char*)tile_out->output_buff;
  size_t output_bytes = graph.output_bytes_per_task;
  std::vector<const char *> input_ptrs;
  std::vector<size_t> input_bytes;
  
  // Add all input dependencies
  for (size_t i = 0; i < tile_in.size(); i++) {
    if (tile_in[i] == NULL || tile_in[i]->output_buff == NULL) {
      fprintf(stderr, "Error: NULL pointer in task2 input %zu at position (%d, %d)\n", 
              i, payload.x, payload.y);
      return;
    }
    input_ptrs.push_back((char*)tile_in[i]->output_buff);
    input_bytes.push_back(graph.output_bytes_per_task);
  }
  
  graph.execute_point(payload.y, payload.x, output_ptr, output_bytes,
                     input_ptrs.data(), input_bytes.data(), input_ptrs.size(), 
                     scratch_memory, graph.scratch_bytes_per_task);
#else  
  tile_out->dep = 0;
  for (size_t i = 0; i < tile_in.size(); i++) {
    if (tile_in[i] == NULL) {
      fprintf(stderr, "Error: NULL pointer in task2 input %zu at position (%d, %d)\n", 
              i, payload.x, payload.y);
      continue;
    }
    tile_out->dep += tile_in[i]->dep;
  }
  printf("Task2 x %d, y %d, out %f\n", payload.x, payload.y, tile_out->dep);
#endif
}

// Similar modifications for task3-task10...

struct SerialApp : public App {
  SerialApp(int argc, char **argv);
  ~SerialApp();
  void execute_main_loop();
  void execute_timestep(size_t idx, long t);
private:
  void debug_printf(int verbose_level, const char *format, ...);
private:
  matrix_t *matrix;
};

SerialApp::SerialApp(int argc, char **argv)
  : App(argc, argv)
{ 
  matrix = (matrix_t *)malloc(sizeof(matrix_t) * graphs.size());
  
  size_t max_scratch_bytes_per_task = 0;
  
  for (unsigned i = 0; i < graphs.size(); i++) {
    TaskGraph &graph = graphs[i];
    
    matrix[i].M = graph.nb_fields;
    matrix[i].N = graph.max_width;
    matrix[i].data = (tile_t*)malloc(sizeof(tile_t) * matrix[i].M * matrix[i].N);
  
    for (int j = 0; j < matrix[i].M * matrix[i].N; j++) {
      matrix[i].data[j].output_buff = (char *)malloc(sizeof(char) * graph.output_bytes_per_task);
    }
    
    if (graph.scratch_bytes_per_task > max_scratch_bytes_per_task) {
      max_scratch_bytes_per_task = graph.scratch_bytes_per_task;
    }
  }
  
  if (max_scratch_bytes_per_task > 0) {
    scratch_memory = (char*)malloc(sizeof(char)*max_scratch_bytes_per_task);
    TaskGraph::prepare_scratch(scratch_memory, sizeof(char)*max_scratch_bytes_per_task);
  }
}

SerialApp::~SerialApp()
{
  for (unsigned i = 0; i < graphs.size(); i++) {
    for (int j = 0; j < matrix[i].M * matrix[i].N; j++) {
      free(matrix[i].data[j].output_buff);
      matrix[i].data[j].output_buff = NULL;
    }
    free(matrix[i].data);
    matrix[i].data = NULL;
  }
  
  free(matrix);
  matrix = NULL;
  
  if (scratch_memory != NULL) {
    free(scratch_memory);
    scratch_memory = NULL;
  }
}

void SerialApp::execute_main_loop()
{ 
  display();
  
  Timer::time_start();
  
  for (unsigned i = 0; i < graphs.size(); i++) {
    const TaskGraph &g = graphs[i];
    for (int y = 0; y < g.timesteps; y++) {
      execute_timestep(i, y);
    }
  }
  
  double elapsed = Timer::time_end();
  report_timing(elapsed);
}

void SerialApp::execute_timestep(size_t idx, long t)
{
  if (idx >= graphs.size()) {
    fprintf(stderr, "Error: Invalid graph index %zu (max: %zu)\n", idx, graphs.size()-1);
    return;
  }

  const TaskGraph &g = graphs[idx];
  long offset = g.offset_at_timestep(t);
  long width = g.width_at_timestep(t);
  long dset = g.dependence_set_at_timestep(t);
  int nb_fields = g.nb_fields;
  
  if (matrix == NULL || matrix[idx].data == NULL) {
    fprintf(stderr, "Error: NULL matrix data for graph %zu\n", idx);
    return;
  }
  
  if (nb_fields <= 0) {
    fprintf(stderr, "Error: Invalid number of fields %d for graph %zu\n", nb_fields, idx);
    return;
  }
  
  payload_t payload;
  payload.graph = g;
  
  for (int x = offset; x <= offset+width-1; x++) {
    printf("timestep: %d, dep: %ld\n", t, x);
    if (x < 0 || x >= matrix[idx].N) {
      fprintf(stderr, "Error: x index %d out of bounds [0, %d) at timestep %ld\n", 
              x, matrix[idx].N, t);
      continue;
    }
    
    // NonKernelSplit();
    std::vector<std::pair<long, long> > deps = g.dependencies(dset, x);   
    
    payload.x = x;
    payload.y = t;
    
    int row_idx = t % nb_fields;
    long array_idx = row_idx * matrix[idx].N + x;
    
    if (array_idx < 0 || array_idx >= matrix[idx].M * matrix[idx].N) {
      fprintf(stderr, "Error: Array index %ld out of bounds [0, %d) at position (%d, %ld)\n", 
              array_idx, matrix[idx].M * matrix[idx].N, x, t);
      continue;
    }
    
    if (deps.size() == 0) {
      // No dependencies, execute task1
      task1(&matrix[idx].data[array_idx], payload);
    } else {
      if (t == 0) {
        // First timestep, execute task1
        task1(&matrix[idx].data[array_idx], payload);
      } else {
        // Handle dependencies
        tile_t *out = &matrix[idx].data[array_idx];
        std::vector<tile_t*> inputs;
        
        for (size_t i = 0; i < deps.size(); i++) {
          long dep_x = deps[i].first;
          
          if (dep_x < 0 || dep_x >= matrix[idx].N) {
            fprintf(stderr, "Error: Dependency x index %ld out of bounds [0, %d) at position (%d, %ld)\n", 
                    dep_x, matrix[idx].N, x, t);
            continue;
          }
          
          int prev_row_idx = (t-1) % nb_fields;
          long dep_array_idx = prev_row_idx * matrix[idx].N + dep_x;
          
          if (dep_array_idx < 0 || dep_array_idx >= matrix[idx].M * matrix[idx].N) {
            fprintf(stderr, "Error: Dependency array index %ld out of bounds [0, %d) at position (%d, %ld)\n", 
                    dep_array_idx, matrix[idx].M * matrix[idx].N, x, t);
            continue;
          }
          
          inputs.push_back(&matrix[idx].data[dep_array_idx]);
        }
        
        task2(out, inputs, payload);
      }
    } 
    // NonKernelSplit();
  }
}

void SerialApp::debug_printf(int verbose_level, const char *format, ...)
{
  if (verbose_level > VERBOSE_LEVEL) {
    return;
  }
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
}

int main(int argc, char ** argv)
{
  argc = 8;
  argv[0] = "-steps";
  argv[1] = "2";
  argv[2] = "-width";
  argv[3] = "2";
  argv[4] = "-type";
  argv[5] = "stencil_1d";
  argv[6] = "-kernel";
  argv[7] = "compute_bound";
  argv[8] = "-iter";
  argv[9] = "4096";
  argv[10] = "-worker";
  argv[11] = "1";
  SerialApp app(argc, argv);
  app.execute_main_loop();
  return 0;
}
}