//! Runtime configuration shared by core engines (filled in later milestones).

/// Process-wide / reader-level runtime options.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RuntimeConfig {
    /// Worker thread hint (`0` means auto / single-thread default until M5).
    pub threads: usize,
    /// Bounded channel capacity multiplier relative to worker count.
    pub queue_capacity_factor: usize,
    /// Whether parallel results must be merged in deterministic order.
    pub ordered: bool,
}

impl Default for RuntimeConfig {
    fn default() -> Self {
        Self {
            threads: 0,
            queue_capacity_factor: 2,
            ordered: true,
        }
    }
}

impl RuntimeConfig {
    /// Suggested bounded queue capacity for producer/consumer channels.
    pub fn queue_capacity(&self, worker_count: usize) -> usize {
        let workers = worker_count.max(1);
        workers.saturating_mul(self.queue_capacity_factor.max(1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_ordered_and_bounded() {
        let cfg = RuntimeConfig::default();
        assert!(cfg.ordered);
        assert_eq!(cfg.queue_capacity(8), 16);
    }
}
