import { CHEM } from "./chem-content.js";
import { mountApp } from "./engine.js";
import { mountExtras } from "./chem-extras.js";
import { mountAccount } from "./account.js";
import { mountStats } from "./stats.js";

mountAccount();
mountApp(CHEM);
mountExtras();
mountStats("chemistry");
