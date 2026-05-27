# Unit tests for `.heartbeat()` — the TTY-aware progress emitter introduced in
# PR #130 follow-up. Both branches are pinned at the byte level so a future
# refactor can't quietly break either contract:
#
#   TTY     (`.tty = TRUE`)  → `\r<msg>\033[K` per emit, `\n` on done = TRUE.
#   Non-TTY (`.tty = FALSE`) → cli `ℹ` message per emit, no-op on done = TRUE.
#
# Tests pass `.tty` explicitly so they don't depend on the test runner's actual
# stderr being (or not being) a terminal. Production call sites always omit
# `.tty` and let `isatty(stderr())` pick the branch — covered by separate
# end-to-end dogfood (not regression-tested in CI).

test_that(".heartbeat TTY mode emits \\r<msg>\\033[K bytes per call", {
  # The canonical in-place pattern: carriage return → message → clear-to-EOL.
  # `\033[K` (ESC[K) clears from cursor to end-of-line so a shorter follow-up
  # message doesn't leave stale chars from a longer prior one. Verified at the
  # byte level — anyone breaking the escape would silently re-introduce the
  # ghosting bug that the `space = paste0(rep(' ',50)...)` workarounds in
  # toposearch.R were originally papering over.
  out = capture.output(.heartbeat("hello", .tty = TRUE), type = "message")
  expect_length(out, 1)
  expect_identical(charToRaw(out[1]),
                   as.raw(c(0x0d, 0x68, 0x65, 0x6c, 0x6c, 0x6f, 0x1b, 0x5b, 0x4b)))
})


test_that(".heartbeat TTY mode glue-interpolates from caller's frame", {
  # The helper grabs parent.frame() so call sites can use the same `{var}`
  # syntax as `cli::cli_inform()`. Tested via a local variable in this frame.
  i = 7
  n = 50
  out = capture.output(.heartbeat("part {i} of {n}", .tty = TRUE),
                       type = "message")
  expect_match(out[1], "part 7 of 50", fixed = TRUE)
  # Still has the \r prefix and \033[K suffix.
  bytes = charToRaw(out[1])
  expect_identical(bytes[1], as.raw(0x0d))
  expect_identical(tail(bytes, 3), as.raw(c(0x1b, 0x5b, 0x4b)))
})


test_that(".heartbeat TTY mode done = TRUE emits a single newline", {
  # The teardown signal that advances the cursor past the heartbeat line.
  out = capture.output(.heartbeat(done = TRUE, .tty = TRUE), type = "message")
  expect_length(out, 1)
  expect_identical(out[1], "")  # capture.output strips the trailing \n,
                                # leaving an empty-string element. The fact
                                # that we got *one* line back (not zero)
                                # confirms a newline was emitted.
})


test_that(".heartbeat TTY mode is silent when msg is empty and not done", {
  # Defensive: `nzchar(msg)` guard prevents an unprovoked `\r\033[K` emit
  # if a caller accidentally passes `""`. We only want to clear the line on
  # an explicit `done = TRUE` teardown.
  out = capture.output(.heartbeat("", .tty = TRUE), type = "message")
  expect_length(out, 0)
})


test_that(".heartbeat TTY mode overwrites in place across successive calls", {
  # Two emits should produce two `\r<msg>\033[K` cycles back-to-back, with
  # no intervening newline. In a real terminal the second emit visually
  # overwrites the first.
  out = capture.output({
    .heartbeat("first",  .tty = TRUE)
    .heartbeat("second", .tty = TRUE)
  }, type = "message")
  bytes = charToRaw(paste(out, collapse = ""))
  # \r first \033[K \r second \033[K
  expected = c(charToRaw("\rfirst"),  as.raw(c(0x1b, 0x5b, 0x4b)),
               charToRaw("\rsecond"), as.raw(c(0x1b, 0x5b, 0x4b)))
  expect_identical(bytes, expected)
})


test_that(".heartbeat non-TTY mode emits a cli ℹ-bulleted message", {
  # The orchestrator-tail contract: each emit lands on stderr as its own
  # newline-terminated line, with cli's ℹ bullet (U+2139) and no ANSI
  # color codes when stderr is not a TTY. cli::cli_inform routes through
  # `rlang::inform()` so testthat's `capture_messages()` (which hooks the
  # condition system) catches it.
  out = testthat::capture_messages(.heartbeat("hello", .tty = FALSE))
  expect_length(out, 1)
  # Bullet glyph + space + msg + \n. The glyph is `ℹ` (U+2139) when cli
  # detects unicode output, and the ASCII fallback `i` otherwise.
  # testthat edition 3 sets `cli.unicode = FALSE` inside each test_that()
  # block (for reproducible output), so under `Config/testthat/edition: 3`
  # we get the ASCII form. Pin the message content rather than the glyph.
  expect_match(out[1], "[ℹi] hello")
})


test_that(".heartbeat non-TTY mode glue-interpolates from caller's frame", {
  i = 3
  n = 10
  out = testthat::capture_messages(
    .heartbeat("part {i} of {n}", .tty = FALSE))
  expect_match(out[1], "part 3 of 10", fixed = TRUE)
})


test_that(".heartbeat non-TTY mode done = TRUE is a silent no-op", {
  # In non-TTY mode each prior emit was already newline-terminated, so the
  # teardown has nothing to do. Tested both via the condition system
  # (capture_messages) and via stderr (capture.output type="message") in
  # case cli ever changed its emission channel.
  msgs = testthat::capture_messages(.heartbeat(done = TRUE, .tty = FALSE))
  expect_length(msgs, 0)
  raw  = capture.output(.heartbeat(done = TRUE, .tty = FALSE), type = "message")
  expect_length(raw, 0)
})


test_that(".heartbeat non-TTY mode is silenced by suppressMessages()", {
  # The cli_inform route emits a `message` condition, which suppressMessages()
  # catches. This is the documented escape hatch for orchestrator runs that
  # want a quieter log. The TTY branch uses cat() to stderr which is NOT a
  # message condition — covered in the next test.
  out = testthat::capture_messages(
    suppressMessages(.heartbeat("noisy", .tty = FALSE)))
  expect_length(out, 0)
})


test_that(".heartbeat TTY mode is NOT silenced by suppressMessages()", {
  # Documented design choice: TTY users who want silence pass verbose = FALSE
  # at the call site. `cat(file = stderr())` is not a message condition, so
  # suppressMessages can't catch it. This test pins that contract — if
  # someone later "fixes" the TTY branch by routing through message(), they'd
  # silently regress the in-place overwrite (message() adds its own newline
  # and prefix, breaking the \r-clear-EOL pattern). The asymmetry is
  # intentional and tested.
  out = capture.output(suppressMessages(.heartbeat("loud", .tty = TRUE)),
                       type = "message")
  expect_length(out, 1)
  expect_match(out[1], "loud", fixed = TRUE)
})
