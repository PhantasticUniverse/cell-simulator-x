//! Screenshot capture functionality.

use std::path::PathBuf;

use anyhow::Result;
use chrono::Local;

/// Save a screenshot to the screenshots directory.
///
/// Creates the screenshots directory if it doesn't exist.
/// Filename is auto-generated with timestamp: `screenshot_YYYYMMDD_HHMMSS.png`
///
/// Returns the path to the saved screenshot.
pub fn save_screenshot(pixels: &[u8], width: u32, height: u32) -> Result<PathBuf> {
    // Create screenshots directory
    let dir = PathBuf::from("screenshots");
    std::fs::create_dir_all(&dir)?;

    // Generate filename with timestamp
    let timestamp = Local::now().format("%Y%m%d_%H%M%S");
    let filename = format!("screenshot_{}.png", timestamp);
    let path = dir.join(&filename);

    // Save image
    // Note: pixels are expected to be in RGBA format
    image::save_buffer(&path, pixels, width, height, image::ColorType::Rgba8)?;

    log::info!("Screenshot saved: {}", path.display());
    Ok(path)
}

/// Save screenshot with custom filename
pub fn save_screenshot_as(pixels: &[u8], width: u32, height: u32, filename: &str) -> Result<PathBuf> {
    // Create screenshots directory
    let dir = PathBuf::from("screenshots");
    std::fs::create_dir_all(&dir)?;

    let path = dir.join(filename);
    image::save_buffer(&path, pixels, width, height, image::ColorType::Rgba8)?;

    log::info!("Screenshot saved: {}", path.display());
    Ok(path)
}
