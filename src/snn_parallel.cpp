#include <Rcpp.h>
#include <RcppParallel.h>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <map>
#include <thread>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

namespace {

int configured_thread_count() {
  const char* configured = std::getenv("RCPP_PARALLEL_NUM_THREADS");
  if (configured != NULL) {
    char* end = NULL;
    const long parsed = std::strtol(configured, &end, 10);
    if (end != configured && *end == '\0' && parsed > 0 &&
        parsed <= std::numeric_limits<int>::max()) {
      return static_cast<int>(parsed);
    }
  }
  const unsigned int available = std::thread::hardware_concurrency();
  return available > 0 ? static_cast<int>(available) : 1;
}

std::size_t worker_grain_size(std::size_t n_cells, int n_threads) {
  const std::size_t threads = static_cast<std::size_t>(n_threads);
  return std::max<std::size_t>(1, (n_cells + threads - 1) / threads);
}

std::vector<int> build_posting_indptr(const RMatrix<int>& knn) {
  const std::size_t n_cells = knn.nrow();
  const std::size_t k = knn.ncol();
  std::vector<int> counts(n_cells, 0);
  for (std::size_t i = 0; i < n_cells; ++i) {
    for (std::size_t t = 0; t < k; ++t) {
      ++counts[static_cast<std::size_t>(knn(i, t))];
    }
  }
  std::vector<int> indptr(n_cells + 1, 0);
  for (std::size_t i = 0; i < n_cells; ++i) {
    indptr[i + 1] = indptr[i] + counts[i];
  }
  return indptr;
}

std::vector<int> build_posting_rows(
    const RMatrix<int>& knn,
    const std::vector<int>& indptr
) {
  const std::size_t n_cells = knn.nrow();
  const std::size_t k = knn.ncol();
  std::vector<int> rows(static_cast<std::size_t>(indptr.back()));
  std::vector<int> next(indptr.begin(), indptr.end() - 1);
  for (std::size_t i = 0; i < n_cells; ++i) {
    for (std::size_t t = 0; t < k; ++t) {
      const int neighbor = knn(i, t);
      rows[static_cast<std::size_t>(next[static_cast<std::size_t>(neighbor)]++)] =
        static_cast<int>(i);
    }
  }
  return rows;
}

std::vector<int> edges_for_row(
    const RMatrix<int>& knn,
    const std::vector<int>& posting_indptr,
    const std::vector<int>& posting_rows,
    std::size_t row,
    int min_shared
) {
  std::map<int, int> shared;
  for (std::size_t t = 0; t < knn.ncol(); ++t) {
    const int neighbor = knn(row, t);
    for (int pos = posting_indptr[static_cast<std::size_t>(neighbor)];
         pos < posting_indptr[static_cast<std::size_t>(neighbor) + 1];
         ++pos) {
      const int other = posting_rows[static_cast<std::size_t>(pos)];
      if (other > static_cast<int>(row)) {
        ++shared[other];
      }
    }
  }
  std::vector<int> out;
  out.reserve(shared.size());
  for (const auto& item : shared) {
    if (item.second >= min_shared) {
      out.push_back(item.first);
    }
  }
  return out;
}

struct CountWorker : public Worker {
  const RMatrix<int> knn;
  const std::vector<int>& posting_indptr;
  const std::vector<int>& posting_rows;
  const int min_shared;
  RVector<int> row_counts;

  CountWorker(const IntegerMatrix knn, const std::vector<int>& posting_indptr,
              const std::vector<int>& posting_rows, int min_shared,
              IntegerVector row_counts)
    : knn(knn), posting_indptr(posting_indptr), posting_rows(posting_rows),
      min_shared(min_shared), row_counts(row_counts) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t row = begin; row < end; ++row) {
      row_counts[row] = static_cast<int>(
        edges_for_row(knn, posting_indptr, posting_rows, row, min_shared).size()
      );
    }
  }
};

struct FillWorker : public Worker {
  const RMatrix<int> knn;
  const std::vector<int>& posting_indptr;
  const std::vector<int>& posting_rows;
  const int min_shared;
  const std::vector<int>& row_offsets;
  RMatrix<int> edges;

  FillWorker(const IntegerMatrix knn, const std::vector<int>& posting_indptr,
             const std::vector<int>& posting_rows, int min_shared,
             const std::vector<int>& row_offsets, IntegerMatrix edges)
    : knn(knn), posting_indptr(posting_indptr), posting_rows(posting_rows),
      min_shared(min_shared), row_offsets(row_offsets), edges(edges) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t row = begin; row < end; ++row) {
      const std::vector<int> retained =
        edges_for_row(knn, posting_indptr, posting_rows, row, min_shared);
      int pos = row_offsets[row];
      for (const int other : retained) {
        edges(pos, 0) = static_cast<int>(row);
        edges(pos, 1) = other;
        ++pos;
      }
    }
  }
};

}  // namespace

// [[Rcpp::export]]
IntegerMatrix bcmp_snn_upper_edges_parallel_cpp(
    const IntegerMatrix nn_idx_k,
    const int min_shared
) {
  if (min_shared < 0) {
    stop("min_shared must be non-negative.");
  }
  const RMatrix<int> knn(nn_idx_k);
  const std::size_t n_cells = knn.nrow();
  if (n_cells == 0 || knn.ncol() == 0) {
    return IntegerMatrix(0, 2);
  }

  const std::vector<int> posting_indptr = build_posting_indptr(knn);
  const std::vector<int> posting_rows = build_posting_rows(knn, posting_indptr);
  IntegerVector row_counts(static_cast<R_xlen_t>(n_cells));
  CountWorker count_worker(
    nn_idx_k, posting_indptr, posting_rows, min_shared, row_counts
  );
  const int n_threads = configured_thread_count();
  const std::size_t grain_size = worker_grain_size(n_cells, n_threads);
  parallelFor(0, n_cells, count_worker, grain_size, n_threads);

  std::vector<int> row_offsets(n_cells + 1, 0);
  for (std::size_t row = 0; row < n_cells; ++row) {
    const int count = row_counts[row];
    if (count > std::numeric_limits<int>::max() - row_offsets[row]) {
      stop("SNN edge count exceeds R integer-matrix capacity.");
    }
    row_offsets[row + 1] = row_offsets[row] + count;
  }

  IntegerMatrix edges(row_offsets.back(), 2);
  FillWorker fill_worker(
    nn_idx_k, posting_indptr, posting_rows, min_shared, row_offsets, edges
  );
  parallelFor(0, n_cells, fill_worker, grain_size, n_threads);
  return edges;
}
