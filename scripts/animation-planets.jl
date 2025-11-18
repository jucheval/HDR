using WGLMakie
using Bonito
using DrWatson

App() do
    # Constants for layout and geometry
    FIGURE_RESOLUTION = (1400, 800)
    TRIANGLE_RADIUS = 2.0             # Distance from the middle box to the outer boxes
    BOX_HEIGHT = 1.0                  # Height of box
    BOX_WIDTH = 16 / 9 * BOX_HEIGHT   # Width of box
    CORNER_RADIUS = 0.1               # Radius of box corners
    ZOOM_PADDING = 1.1                # Padding factor around zoomed box
    ANIMATION_DURATION = 2.0          # Animation duration in seconds
    FRAMES_PER_SECOND = 60            # Animation frames per second
    EASING_FUNCTION = t -> 3t^2 - 2t^3  # Smooth step easing (cubic Hermite)
    BOX_LABELS = ("""
        Point processes
        Neurosciences
        """,
        "Probability",
        "PDE",
        "Statistics")

    # Function to draw a rounded rectangle box
    function RoundedRectangle(center::Point2, width::Real, height::Real, radius::Real)
        # Clamp radius to valid range
        r = clamp(Float64(radius), 0.0, min(width, height) / 2)
        # retrieve coordinates
        xc, yc = center
        # compute bottomleft coordinates
        x, y = (xc - width / 2, yc - height / 2)

        # Zero radius? Use a standard rectangle
        if r < 1e-10
            return Rect(x, y, width, height)
        end

        # Bezier control point offset for quarter circle approximation
        # This magic number is 4/3 * (sqrt(2) - 1) ≈ 0.5522847498
        k = 0.5522847498 * r

        # Build the path (clockwise from top edge)
        return BezierPath([
            MoveTo(Point2f(x + width - r, y)),
            LineTo(Point2f(x + r, y)),
            CurveTo(Point2f(x + r - k, y), Point2f(x, y + r - k), Point2f(x, y + r)),
            LineTo(Point2f(x, y + height - r)),
            CurveTo(Point2f(x, y + height - r + k), Point2f(x + r - k, y + height), Point2f(x + r, y + height)),
            LineTo(Point2f(x + width - r, y + height)),
            CurveTo(Point2f(x + width - r + k, y + height), Point2f(x + width, y + height - r + k), Point2f(x + width, y + height - r)),
            LineTo(Point2f(x + width, y + r)),
            CurveTo(Point2f(x + width, y + r - k), Point2f(x + width - r + k, y), Point2f(x + width - r, y)),
            ClosePath()
        ])
    end

    # Animation function with lock to prevent concurrent animations
    function animate_limits!(target_x, target_y; zoom_in=true)
        global is_animating, current_xlims, current_ylims

        if is_animating
            return
        end

        # Does nothing if already at target
        if target_x == current_xlims && target_y == current_ylims

            return
        end

        is_animating = true

        # Calculate animation parameters
        total_frames = ceil(Int, ANIMATION_DURATION * FRAMES_PER_SECOND)
        middle_frame = div(total_frames, 2)
        start_xlims = current_xlims
        start_ylims = current_ylims

        for frame in 1:total_frames
            # Check if figure is still open
            if !isopen(fig.scene)
                is_animating = false
                return
            end

            t_raw = frame / total_frames
            t_switch_in = 0.6
            t_switch_out = 1.0 - t_switch_in

            if zoom_in
                alpha_background[] = EASING_FUNCTION(t_raw / t_switch_in) * (t_raw <= t_switch_in) + EASING_FUNCTION(1.0 - (t_raw - t_switch_in) / (1.0 - t_switch_in)) * (t_raw > t_switch_in)
                notify(alpha_background)
                t_raw == t_switch_in && (visibles_labels[] = false; notify(visibles_labels))
                t_raw == t_switch_in && (visible_images[] = true; notify(visible_images))
            end
            if !zoom_in
                alpha_background[] = EASING_FUNCTION(t_raw / t_switch_out) * (t_raw <= t_switch_out) + EASING_FUNCTION(1.0 - (t_raw - t_switch_out) / (1.0 - t_switch_out)) * (t_raw > t_switch_out)
                notify(alpha_background)
                t_raw == t_switch_out && (visibles_labels[] = true; notify(visibles_labels))
                t_raw == t_switch_out && (visible_images[] = false; notify(visible_images))
            end

            # Apply easing function for smooth acceleration/deceleration
            t_smooth = EASING_FUNCTION(t_raw)

            # Interpolate limits
            new_xlims = (
                start_xlims[1] + (target_x[1] - start_xlims[1]) * t_smooth,
                start_xlims[2] + (target_x[2] - start_xlims[2]) * t_smooth
            )
            new_ylims = (
                start_ylims[1] + (target_y[1] - start_ylims[1]) * t_smooth,
                start_ylims[2] + (target_y[2] - start_ylims[2]) * t_smooth
            )

            # Update axis limits
            xlims!(ax, new_xlims)
            ylims!(ax, new_ylims)

            # Update current state
            current_xlims = new_xlims
            current_ylims = new_ylims

            # Frame timing
            sleep(1 / FRAMES_PER_SECOND)
        end

        # Ensure final position is exact
        xlims!(ax, target_x)
        ylims!(ax, target_y)
        current_xlims = target_x
        current_ylims = target_y

        is_animating = false
    end

    function animate_limits_2steps!(target_x, target_y)
        global is_animating, current_xlims, current_ylims

        # # Prevent concurrent animations
        # if is_animating
        #     return
        # end
        # #is_animating = true

        @async begin
            # Reset view
            animate_limits!(initial_xlims, initial_ylims; zoom_in=false)
            # Move view
            animate_limits!(target_x, target_y; zoom_in=true)
        end
    end

    begin # Initialize figure
        # Create figure and axis
        fig = Figure(size=FIGURE_RESOLUTION)
        ax = Axis(fig[1, 2], aspect=DataAspect())
        hidedecorations!(ax)
        hidespines!(ax)

        # Centers of circles
        centers = [
            Point2f(0, 0),
            Point2f(0, TRIANGLE_RADIUS),
            Point2f(-TRIANGLE_RADIUS, -0.7 * TRIANGLE_RADIUS),
            Point2f(TRIANGLE_RADIUS, -0.7 * TRIANGLE_RADIUS)
        ]

        # Strokecolors and labels for the circles
        colors = [:blue, :green, :red, :purple]
        visibles_labels = Observable(true)
        visible_images = Observable(false)
        alpha_background = Observable(0.0)

        for i in 1:4
            # Box content as image
            xpad = BOX_WIDTH / 2
            ypad = BOX_HEIGHT / 2
            img_x = centers[i][1] .+ (-1.0, +1.0) .* xpad
            img_y = centers[i][2] .+ (-1.0, +1.0) .* ypad
            img = load("notebooks/assets/box" * string(i) * ".png")
            image!(ax, img_x, img_y, rotr90(img), visible=visible_images)

            # Box label
            text!(ax, centers[i],
                text=BOX_LABELS[i],
                align=(:center, :center),
                fontsize=28,
                color=:black,
                font=:bold,
                visible=visibles_labels)

            # Box background to deal with transparency
            poly!(ax, RoundedRectangle(centers[i], BOX_WIDTH, BOX_HEIGHT, CORNER_RADIUS),
                color=:white,
                alpha=alpha_background)

            # Box contour
            poly!(ax, RoundedRectangle(centers[i], BOX_WIDTH, BOX_HEIGHT, CORNER_RADIUS),
                color=:transparent,
                strokecolor=Makie.Colors.JULIA_LOGO_COLORS[colors[i]],
                strokewidth=4)
        end

        # Set initial limits to show all boxes
        x_extent = TRIANGLE_RADIUS + BOX_WIDTH / 2
        y_extent = TRIANGLE_RADIUS + BOX_HEIGHT / 2
        initial_xlims = (-x_extent * 1.1, x_extent * 1.1)
        initial_ylims = (-y_extent * 1.1, y_extent * 1.1)
        xlims!(ax, initial_xlims)
        ylims!(ax, initial_ylims)

        # Create selection buttons
        buttongrid = fig[1, 1] = GridLayout(tellheight=false)
        buttonlabels = ["Show all", "PP/Neuro", "Probability", "PDE", "Statistics"]
        buttons = buttongrid[1:5, 1] = [Makie.Button(fig, label=buttonlabels[i], width=Fixed(100)) for i in 1:5]

        # Animation state
        global is_animating = false
        global current_xlims = initial_xlims
        global current_ylims = initial_ylims
    end

    # Connect buttons to zoom action
    for buttonid in 1:5
        on(buttons[buttonid].clicks) do _
            if buttonid == 1
                target_x = initial_xlims
                target_y = initial_ylims
            else
                circle_num = buttonid - 1
                center = centers[circle_num]

                xpad = BOX_WIDTH * ZOOM_PADDING / 2
                ypad = BOX_HEIGHT * ZOOM_PADDING / 2
                target_x = (center[1] - xpad, center[1] + xpad)
                target_y = (center[2] - ypad, center[2] + ypad)
            end

            # Start animation
            animate_limits_2steps!(target_x, target_y)
        end
    end

    # Adjust layout for better appearance
    colsize!(fig.layout, 1, Fixed(150))

    DOM.div(fig)
end





