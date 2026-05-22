use std::fs;
use std::fs::OpenOptions;
use std::io::Write;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

pub fn unique_temp_dir(prefix: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("{}_{}", prefix, nanos));
    fs::create_dir_all(&dir).expect("failed to create temp dir");
    dir
}

pub fn repo_log_path(prefix: &str) -> PathBuf {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let log_dir = manifest_dir.join("test_output/integration_logs");
    fs::create_dir_all(&log_dir).expect("failed to create integration log dir");
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before unix epoch")
        .as_nanos();
    log_dir.join(format!("{}_{}.log", prefix, nanos))
}

pub fn log_line(log_path: &std::path::Path, message: &str) {
    let mut log = OpenOptions::new()
        .create(true)
        .append(true)
        .open(log_path)
        .expect("failed to open integration log");
    writeln!(log, "{}", message).expect("failed to write integration log");
}

pub fn log_command(log_path: &std::path::Path, binary: &str, args: &[&str]) {
    let cmd_str = format!("{} {}", binary, args.join(" "));
    log_line(log_path, &format!("Executing: {}", cmd_str));
}
