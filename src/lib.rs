mod data;
mod environment;
mod prelude;
mod solver;
mod state;
mod trajectory;

use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    fn alert(s: &str);
}