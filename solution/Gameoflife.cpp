#include <algorithm>
#include <cassert>
#include <exception>
#include <iostream>
#include <vector>

#include "GUI.h"
#include "Gameoflife.h"
#include "Graph.h"
#include "Window.h"

#include "FL/fl_ask.H"

using namespace Graph_lib;

Cell::Cell(Point pos, int size)
	: state(State::Dead), rect{std::make_shared<Rectangle>(pos, size, size)}
{
	// Do not show the border of the cell
	rect->set_color(border_hide);

	// Update the graphical state based on the current state
	update();
}

// TASK C1
int Cell::get_value() const
{
	// BEGIN: C1
	return static_cast<int>(state);
	// BEGIN: C1
}

// TASK C2
void Cell::update()
{
	// BEGIN: C2
	rect->set_fill_color(colors[get_value()]);
	// BEGIN: C2
}

// TASK C3
void Cell::kill()
{
	// BEGIN: C3
	state = State::Dead;
	update();
	// END: C3
}

// TASK C4
void Cell::resurrect()
{
	// BEGIN: C4
	state = State::Live;
	update();
	// END: C4
}

// TASK C5
void Cell::set_state(char c)
{
	// BEGIN: C5
	state = char_to_state.at(c);
	update();
	// END: C5
}

// TASK C6
std::istream& operator>>(std::istream& is, Cell& cell)
{
	// BEGIN: C6
	char state;
	is >> state;
	cell.set_state(state);
	return is;
	// END: C6
}

// TASK C7
bool Cell::is_alive() const
{
	// BEGIN: C7
	return state == State::Live;
	// END: C7
}

// TASK C8
char Cell::as_char() const
{
	// BEGIN: C8
	return chars[get_value()];
	// BEGIN: C8
}

// TASK G1
Gameoflife::Gameoflife(int x_cells,
                       int y_cells,
                       const std::string& window_title)
	: Window{x_cells * cell_size + 4 * button_width + 2 * margin,
             y_cells * cell_size + 2 * margin,
             window_title},
	  x_cells{x_cells},
	  y_cells{y_cells}
{
	// Provided (asserts and attach() calls)
	// Asserts are here to make sure any unintended changes to
	// compile-time constants makes the exercise less hard to debug.
	assert(x_cells > 0);
	assert(y_cells > 0);
	assert(button_width > 0);
	assert(margin > 0);

	attach(step_btn);
	attach(steps_input);

	attach(load_btn);
	attach(save_btn);
	attach(filename_input);

	attach(toggle_gridlines_btn);
	attach(play_pause_btn);

	// Provided
	// We do not want the window or its children widgets (button, inputs,
	// etc.) to be resizable. This only leads to possible confusion and
	// frustration. If you struggle with a window that is too large/too
	// small, modify the `cell_size` variable at the top of Gameoflife.h so
	// that it suits your screen resolution and size.
	resizable(nullptr);
	size_range(x_max(), y_max(), x_max(), y_max());

	// BEGIN: G1
	// Populate grids.
	// We use emplace_back() in the solution guide.
	// Solutions with push_back() are equally good for this exam.
	for (int y = 0; y < y_cells; ++y) {
		get_current_grid().emplace_back();
		for (int x = 0; x < x_cells; ++x) {
			int x_pos = x * cell_size + margin;
			int y_pos = y * cell_size + margin;
			get_current_grid()
				.back()
				.emplace_back(Point{x_pos, y_pos}, cell_size)
				.attach_to(*this);
		}
	}

	// Make two identical copies of the grid
	get_scratch_grid() = get_current_grid();

	// For convenience, initially load the flower
	load("flower.cgol");
	// END: G1
}

// TASK G2
std::istream& operator>>(std::istream& is, Gameoflife& gameoflife)
{
	// BEGIN: G2
	for (auto& row : gameoflife.get_current_grid()) {
		for (auto&& cell : row) {
			is >> cell;
		}
	}

	return is;
	// END: G2
}

// TASK G3
void Gameoflife::load(const std::string& filename)
{
	// BEGIN: G3
	std::ifstream ifs{filename};
	if (!ifs) {
		throw std::runtime_error{"Could not load a Game of life state from '" +
		                         filename + "'."};
	}

	ifs >> *this;

	redraw();
	// END: G3
}

// TASK G4
void Gameoflife::step()
{
	// BEGIN: G4
	for (size_t row = 0; row < get_current_grid().size(); ++row) {
		for (size_t col = 0; col < get_current_grid().back().size(); ++col) {

			// How many neighbours are currently alive?
			int live_neighbours = 0;
			// We go through all neighbouring cells
			for (int y = -1; y <= 1; ++y) {
				for (int x = -1; x <= 1; ++x) {
					// Except the current cell
					if (!(x == 0 && y == 0)) {
						// And add the integer value of the cell's state
						int lookup_y = ((row + y) + y_cells) % y_cells;
						int lookup_x = ((col + x) + x_cells) % x_cells;
						live_neighbours +=
							get_current_grid()[lookup_y][lookup_x].get_value();
					}
				}
			}

			// Alternatives
			// Helper lambda to retrieve the value of a cell.
			// This makes the lookup code for the following alternatives
			// easier to read for sensors.
			// auto grid_value = [row, col, this](int x, int y) {
			// We wrap on edges
			// 	int lookup_y = ((row + y) + y_cells) % y_cells;
			// 	int lookup_x = ((col + x) + x_cells) % x_cells;
			// 	return get_current_grid()[lookup_y][lookup_x].get_value();
			// };

			// Kernel computation - How many of the neighbours are currently
			// alive?
			// int live_neighbours =
			// 	grid_value(-1, -1) + grid_value(0, -1) + grid_value(1, -1) +
			// 	grid_value(-1, 0) /*+ grid_value(0,0)*/ + grid_value(1, 0) +
			// 	grid_value(-1, 1) + grid_value(0, 1) + grid_value(1, 1);

			// Alternative with lookup and loop:
			// int live_neighbours = 0;
			// std::array<int, 3> kernel{-1, 0, 1};
			// for (int y : kernel) {
			// 	for (int x : kernel) {
			// 		if (!(x == 0 && y == 0)) {
			// 			live_neighbours += grid_value(x, y);
			// 		}
			// 	}
			// }

			// Yet another one is to define a kernel, then loop through
			// similar to above, unroll, etc.

			// The rules
			if (get_current_grid()[row][col].is_alive() &&
			    (live_neighbours == 2 || live_neighbours == 3)) {
				// 1. Any live cell with two or three live neighbours
				// survives.
				get_scratch_grid()[row][col].resurrect();
				// Pass - keep alive
			} else if (!get_current_grid()[row][col].is_alive() &&
			           live_neighbours == 3) {
				// 2. Any dead cell with three live neighbours becomes a
				// live cell.
				get_scratch_grid()[row][col].resurrect();
			} else {
				// 3. All other live cells die in the next generation.
				// Similarly, all other dead cells stay dead.
				get_scratch_grid()[row][col].kill();
			}
		}
	}

	// Swap the state grids to reflect the newly updated state.
	std::swap(current_grid, scratch_grid);
	// END: G4
}

// TASK G5
void Gameoflife::step(int steps)
{
	// BEGIN: G5
	for (; steps--;) {
		step();
	}

	redraw();
	// END: G5
}

// TASK E1
Cell* Gameoflife::cell_at_pos(Point pos)
{
	// BEGIN: E1
	int row = (pos.y - margin) / cell_size;
	int col = (pos.x - margin) / cell_size;
	std::cout << row << ", " << col << '\n';

	// If the position is outside the grid, we
	// have no cell at the position.
	if (row >= y_cells || col >= x_cells || pos.y - margin < 0 ||
	    pos.x - margin < 0)
		return nullptr;

	return &get_current_grid()[row][col];
	// END: E1
}

// TASK E2
void Cell::toggle()
{
	// BEGIN: E2
	if (state == State::Live) {
		kill();
	} else {
		resurrect();
	}
	// END: E2
}

// TASK E3
bool Gameoflife::toggle_cell(Point pos)
{
	// BEGIN: E3
	if (auto cell = cell_at_pos(pos)) {
		cell->toggle();
		redraw();
		return true;
	}

	return false;
	// END: E3
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
// The following code is provided as part of the exam handout
// and should _not_ be modified.

// Provided
void Cell::toggle_border()
{
	// Either visible or invisible border lines,
	// based on current visibility state
	Color new_grid_color = show_border ? border_hide : border_show;
	rect->set_color(new_grid_color);
	show_border = !show_border;
}

// Provided - same as ex7
void Cell::attach_to(Window& win) { win.attach(*rect); }

// Provided
std::ostream& operator<<(std::ostream& os, const Cell& cell)
{
	return os << cell.as_char();
}

// Provided
void Gameoflife::save(const std::string& filename)
{
	// We do not want anyone to overwrite their
	// source code files. Make sure the last characters
	// are ".cgol".
	constexpr std::string_view ext = ".cgol";
	if (filename.size() < ext.size() ||
	    filename.substr(filename.size() - 5, 5) != ext) {
		throw std::runtime_error{"'" + filename +
		                         "' does not have the correct file extension: "
		                         "'.cgol' or is too short."};
	}

	std::ofstream ofs{filename};
	if (!ofs) {
		throw std::runtime_error{"Could not save the CGoL state to '" +
		                         filename + "'."};
	}

	ofs << *this;
}

// Provided
std::ostream& operator<<(std::ostream& os, const Gameoflife& gameoflife)
{
	for (auto& row : gameoflife.get_current_grid()) {
		for (auto&& cell : row) {
			os << cell;
		}
		os << '\n';
	}

	return os;
}

// Provided
void Gameoflife::step_pressed()
{
	int steps = std::clamp(steps_input.get_int(), min_steps, max_steps);
	step(steps);
}

// Provided
void Gameoflife::load_state()
{
	std::string filename = filename_input.get_string();
	try {
		load(filename);
	} catch (const std::runtime_error& e) {
		alert(e.what());
	}
}

// Provided
void Gameoflife::save_state()
{
	std::string filename = filename_input.get_string();
	try {
		save(filename);
	} catch (const std::runtime_error& e) {
		alert(e.what());
	}
}

// Provided
Gameoflife::Grid& Gameoflife::get_current_grid() { return grid[current_grid]; }
Gameoflife::Grid& Gameoflife::get_scratch_grid() { return grid[scratch_grid]; }
const Gameoflife::Grid& Gameoflife::get_current_grid() const
{
	return grid[current_grid];
}
const Gameoflife::Grid& Gameoflife::get_scratch_grid() const
{
	return grid[scratch_grid];
}

// Provided
// Show a pop-up alert box with `message`
void Gameoflife::alert(const std::string& message)
{
	fl_alert("%s", message.c_str());
}

// Provided
// Show or hide the grid lines
void Gameoflife::toggle_gridlines()
{
	for (auto& row : get_current_grid()) {
		for (auto& cell : row) {
			cell.toggle_border();
		}
	}

	redraw();
}

// Provided
// Handle events.
// If it's a mouse button click event outside Widgets,
// try to handle the event with grid_click()
int Gameoflife::handle(int event)
{
	if (int handled = Fl_Group::handle(event) > 0) {
		// Event has been passed on to a child widget
		return handled;
	} else if (event == FL_PUSH) {
		// Handle mouse button clicks outside the widgets
		return click_handler(Point{Fl::event_x(), Fl::event_y()});
	} else {
		return Fl_Widget::handle(event);
	}
}

// Provided
// The Gameoflife instance's mouse-click handler.
int Gameoflife::click_handler(Point pos) { return toggle_cell(pos); }

// Provided
void Gameoflife::play_pause()
{
	if (playing.load()) {
		playing.store(false);
		worker.join();
	} else {
		playing.store(true);
		worker = thread(&Gameoflife::animate, ref(*this));
	}
}

// Provided
void Gameoflife::animate()
{
	while (playing.load()) {
		Fl::lock();
		step(1);
		Fl::unlock();
		Fl::awake();

		this_thread::sleep_for(animation_interval);
	}
}

// Provided
Gameoflife::~Gameoflife()
{
	// Clean up animation thread if it is playing to avoid
	// ABI abort trap messages.
	playing.store(false);
	if (worker.joinable())
		worker.join();
}
