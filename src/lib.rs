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

#[wasm_bindgen]
pub fn greet() {
    alert("Hello, {{project-name}}!");
}

#[wasm_bindgen]
pub fn array_test(arr: &mut [f64]) {
    for x in arr.iter_mut() {
        *x = 1.0;
    }
}