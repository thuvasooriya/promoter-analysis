#let blue = rgb("#648fff")
#let purple = rgb("#555ef0").darken(50%)
#let magenta = rgb("#dc267f")
#let brown = rgb("#fe6100").darken(50%)
#let yellow = rgb("#ffb000")
#let orange = none
#let red = none
#let green = none

#let individual_assignment_title(
  title,
  author,
  index,
  subtitle,
  course_id,
  course_name,
  department,
  faculty,
  university,
  body,
) = {
  set document(title: title, author: author)
  set page(
    paper: "a4",
    // margin: (x: 1.8cm, y: 1.5cm),
  )
  set text(
    // font: "Times New Roman",
    // size: 12pt
  )
  set par(
    justify: true,
  )

  align(center + top, image("./assets/uom.png", width: 30%))
  align(center, text(17pt)[
    *#title*])
  align(center, text(13pt)[_ #subtitle _])
  align(center, text(14pt)[*#author #index*])
  align(center, text(12pt)[#(
    datetime.today().display("[day] [month repr:long] [year]")
  )])
  align(center + bottom, image("./assets/entc.png", width: 20%))
  align(center + bottom, text(
    11pt,
  )[Submitted in partial fulfillment of the requirements for the module *#course_id* : *#course_name* from _ #department, #faculty, #university _])

  pagebreak(weak: false)
  body
}
