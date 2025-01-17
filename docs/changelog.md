# Changelog

<!-- prettier-ignore -->
Name | PR | Developer | Date | Version
--- | --- | --- | --- | ---
⬇️ Upper bound scipy for compatibility with fcsparser | [38](https://github.com/laminlabs/readfcs/pull/38) | [Koncopd](https://github.com/Koncopd) | 2024-10-21 | 1.1.9
escape $ properly | [37](https://github.com/laminlabs/readfcs/pull/37) | [Zethson](https://github.com/Zethson) | 2024-05-02 |
🏷️ Save object as category | [36](https://github.com/laminlabs/readfcs/pull/36) | [sunnyosun](https://github.com/sunnyosun) | 2024-04-19 | 1.1.8
🍱 Added oetjen18 files | [35](https://github.com/laminlabs/readfcs/pull/35) | [sunnyosun](https://github.com/sunnyosun) | 2023-10-03 | 1.1.7
👷 Remove nox cache | [34](https://github.com/laminlabs/readfcs/pull/34) | [sunnyosun](https://github.com/sunnyosun) | 2023-08-25 |
👷 Test py3.11 | [33](https://github.com/laminlabs/readfcs/pull/33) | [sunnyosun](https://github.com/sunnyosun) | 2023-08-25 | 1.1.6
🚑️ Raise error if the spill index doesn't match channels | [32](https://github.com/laminlabs/readfcs/pull/32) | [sunnyosun](https://github.com/sunnyosun) | 2023-08-24 |
✏️ Fix blank stripping | [30](https://github.com/laminlabs/readfcs/pull/30) | [sunnyosun](https://github.com/sunnyosun) | 2023-08-08 | 1.1.5
🚑️ Make sure channels are sorted by n | [29](https://github.com/laminlabs/readfcs/pull/29) | [sunnyosun](https://github.com/sunnyosun) | 2023-07-07 | 1.1.4
🩹 Remove whitespace only markers | [26](https://github.com/laminlabs/readfcs/pull/26) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-27 | 1.1.3
➖ Remove nbproject from dependencies | [25](https://github.com/laminlabs/readfcs/pull/25) | [Koncopd](https://github.com/Koncopd) | 2023-06-04 | 1.1.2
💚 Fix docs | [23](https://github.com/laminlabs/readfcs/pull/23) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-07 |
🩹 Fill nan with empty string in channels | [22](https://github.com/laminlabs/readfcs/pull/22) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-07 | 1.1.1
🚑 Wrote new channels formater | [21](https://github.com/laminlabs/readfcs/pull/21) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-20 | 1.1.0
👷 Extend CI to py3.8-3.10 | [20](https://github.com/laminlabs/readfcs/pull/20) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-12 | 1.0.4
🩹 Convert list columns to str so it saves properly | [19](https://github.com/laminlabs/readfcs/pull/19) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-08 | 1.0.3
🩹 Make sure `anndata` index is str | [18](https://github.com/laminlabs/readfcs/pull/18) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-06 | 1.0.2
🐛 Change `_channel_names_` to list for saving error | [17](https://github.com/laminlabs/readfcs/pull/17) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-30 | 1.0.1
🎨 Get coverage to 100% | [16](https://github.com/laminlabs/readfcs/pull/16) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-22 | 1.0.0
♻️ Rewrote readfcs based on fcsparser | [15](https://github.com/laminlabs/readfcs/pull/15) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-22 | 0.1.7
🍱 Added another dataset | [14](https://github.com/laminlabs/readfcs/pull/14) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-22 | 0.1.6
🎨 Improve coverage | [13](https://github.com/laminlabs/readfcs/pull/13) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-22 |
♻️ Clean up | [12](https://github.com/laminlabs/readfcs/pull/12) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-22 |
👷 Add codecov | [11](https://github.com/laminlabs/readfcs/pull/11) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-19 |
📝 Make the tutorials executed | [10](https://github.com/laminlabs/readfcs/pull/10) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-19 |
🩹 Deals with spill is None | [9](https://github.com/laminlabs/readfcs/pull/9) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-09 | 0.1.5
🎨 Set default to redindex | [8](https://github.com/laminlabs/readfcs/pull/8) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-07 | 0.1.4
🎨 Fixed index mismatch | [7](https://github.com/laminlabs/readfcs/pull/7) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-07 | 0.1.3
✨ Added a read function | [6](https://github.com/laminlabs/readfcs/pull/6) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-01 | 0.1.2
🔖 Setup changelog and release v0.1.1 | [5](https://github.com/laminlabs/readfcs/pull/5) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-25 |
✨ Allow to convert `FCSFile` to `AnnData` | [1380806](https://github.com/laminlabs/readfcs/commit/f883805feec636a5b160fe57f7b279fa292652a2) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-11 | 0.1.1
✨ Added write_fcs | [3](https://github.com/laminlabs/readfcs/pull/3) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-11 | 0.1.0
