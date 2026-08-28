# Cache Benchmark Results

Payload sizes: 64 KiB, 128 KiB, and 256 KiB. Results are mean bandwidth over three repetitions.

## Filesystem versus SQLite

| Workload | Payload | Filesystem | SQLite | Winner |
|---|---:|---:|---:|---|
| Durable write | 64 KiB | 93.8 MiB/s | 105.5 MiB/s | SQLite 1.12x |
| Durable write | 128 KiB | 108.8 MiB/s | 103.6 MiB/s | Filesystem 1.05x |
| Durable write | 256 KiB | 124.0 MiB/s | 104.4 MiB/s | Filesystem 1.19x |
| Warm read | 64 KiB | 1573.9 MiB/s | 2013.6 MiB/s | SQLite 1.28x |
| Warm read | 128 KiB | 1757.1 MiB/s | 1859.6 MiB/s | SQLite 1.06x |
| Warm read | 256 KiB | 1691.8 MiB/s | 1811.1 MiB/s | SQLite 1.07x |
| 95/5 mixed | 64 KiB | 1286.8 MiB/s | 1562.4 MiB/s | SQLite 1.21x |
| 95/5 mixed | 128 KiB | 1527.2 MiB/s | 1595.8 MiB/s | SQLite 1.04x |
| 95/5 mixed | 256 KiB | 1584.1 MiB/s | 1550.8 MiB/s | Effectively tied |

Submit-only writes favored the filesystem by 56-95x, but those measurements intentionally exclude queue drain and durability. They measure API return latency rather than storage throughput.

## Separate Databases versus One Shared Database

Here, **shared** means `separate_bins=false`: multiple logical cache bins share one SQLite database.

| Workload | Payload | 2 bins | 8 bins | 16 bins |
|---|---:|---:|---:|---:|
| Durable write | 64 KiB | Separate 1.16x | Shared 1.18x | Shared 1.65x |
| Durable write | 128 KiB | Separate 1.28x | Shared 1.18x | Shared 1.50x |
| Durable write | 256 KiB | Separate 1.20x | Shared 1.37x | Shared 1.55x |
| Warm read | 64 KiB | Separate 1.14x | Separate 1.08x | Separate 1.23x |
| Warm read | 128 KiB | Tie | Separate 1.13x | Separate 1.06x |
| Warm read | 256 KiB | Tie | Separate 1.10x | Separate 1.03x |
| 95/5 mixed | 64 KiB | Tie | Shared 1.12x | Shared 1.26x |
| 95/5 mixed | 128 KiB | Shared 1.12x | Shared 1.08x | Shared 1.37x |
| 95/5 mixed | 256 KiB | Tie | Shared 1.06x | Shared 1.33x |

Near-1.0x results should be treated as ties. Most benchmark results had a coefficient of variation below 10%.
