//! Shared base-bucket lookup for depth / pileup hot loops.
//!
//! Buckets: 0=A 1=C 2=G 3=T 4=N(other). 255 = skip (not a base-carrying op).

/// ASCII → bucket lookup table (depth path).
pub(crate) const BASE_BUCKET: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};
