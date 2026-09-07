import { CHEM } from "./chem-content.js";
import { mountApp } from "./engine.js";
import { mountExtras } from "./chem-extras.js";
import { mountAccount } from "./account.js";
import { mountStats } from "./stats.js";
import { mountHome } from "./home.js";
import { mountStudy } from "./study.js";
import { prefs } from "./prefs.js";

mountAccount();
const study = mountStudy(CHEM);
const ui = mountApp(CHEM, study);
mountExtras();
mountStats("chemistry");
mountHome(CHEM, ui, study);
prefs.sync();
