/* Copyright 2020 Los Alamos National Laboratory
 * Licensed under the Apache License, Version 2.0 */

#include <stdarg.h>
#include <assert.h>
#include <string.h>
#include <algorithm> 
#include <unistd.h>
#include "core.h"
#include "timer.h"
#include <map>
#define VERBOSE_LEVEL 0
#define USE_CORE_VERIFICATION

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
#if defined (USE_CORE_VERIFICATION)    
  TaskGraph graph = payload.graph;
  char *output_ptr = (char*)tile_out->output_buff;
  size_t output_bytes = graph.output_bytes_per_task;
  std::vector<const char *> input_ptrs;
  std::vector<size_t> input_bytes;
  
  // Add all input dependencies
  for (size_t i = 0; i < tile_in.size(); i++) {
    input_ptrs.push_back((char*)tile_in[i]->output_buff);
    input_bytes.push_back(graph.output_bytes_per_task);
  }
  
  graph.execute_point(payload.y, payload.x, output_ptr, output_bytes,
                     input_ptrs.data(), input_bytes.data(), input_ptrs.size(), 
                     scratch_memory, graph.scratch_bytes_per_task);
#else  
  tile_out->dep = 0;
  for (size_t i = 0; i < tile_in.size(); i++) {
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
  const TaskGraph &g = graphs[idx];
  long offset = g.offset_at_timestep(t);
  long width = g.width_at_timestep(t);
  long dset = g.dependence_set_at_timestep(t);
  int nb_fields = g.nb_fields;
  
  payload_t payload;
  payload.graph = g;
  
  for (int x = offset; x <= offset+width-1; x++) {
    std::vector<std::pair<long, long> > deps = g.dependencies(dset, x);   
    
    payload.x = x;
    payload.y = t;
    
    if (deps.size() == 0) {
      // No dependencies, execute task1
      task1(&matrix[idx].data[t % nb_fields * matrix[idx].N + x], payload);
    } else {
      if (t == 0) {
        // First timestep, execute task1
        task1(&matrix[idx].data[t % nb_fields * matrix[idx].N + x], payload);
      } else {
        // Handle dependencies
        tile_t *out = &matrix[idx].data[t % nb_fields * matrix[idx].N + x];
        std::vector<tile_t*> inputs;
        for (size_t i = 0; i < deps.size(); i++) {
          inputs.push_back(&matrix[idx].data[(t-1) % nb_fields * matrix[idx].N + deps[i].first]);
        }
        task2(out, inputs, payload);
      }
    }
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
  std::map<uint64_t, uint64_t>* params =  new std::map<uint64_t, uint64_t>();
  params->insert(std::make_pair(1, 10));
  params->insert(std::make_pair(2, 20));
  params->insert(std::make_pair(3, 30));
  params->insert(std::make_pair(4, 40));
  params->insert(std::make_pair(5, 50));
  params->insert(std::make_pair(6, 60));
  params->insert(std::make_pair(7, 70));
  SerialApp app(argc, argv);
  app.execute_main_loop();
  delete params;
  return 0;
}