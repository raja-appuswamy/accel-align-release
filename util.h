#pragma once

static inline int compute_elapsed(struct timeval *start, struct timeval *end) {
  return ((end->tv_sec * 1000000 + end->tv_usec)
      - (start->tv_sec * 1000000 + start->tv_usec));
}